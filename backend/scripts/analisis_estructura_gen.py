#!/usr/bin/env python3
"""
GenomeHub - Analisis de Estructura Genomica Bacteriana

Analiza la composicion del genoma:
- Regiones codificantes vs no codificantes
- Identificacion de operones putativos
- Datos educativos sobre estructura del gen bacteriano
"""

import os
import sys
import json
import gc
from datetime import datetime
from collections import Counter

try:
    from Bio import SeqIO
    from Bio.SeqUtils import gc_fraction
except ImportError as e:
    print(f"[ERROR] BioPython no instalado: {e}")
    sys.exit(1)

# =============================================================================
# CONFIGURACION - CON VALIDACION
# =============================================================================

try:
    DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # backend/
    DIRECTORIO_PROYECTO = os.path.dirname(DIRECTORIO_BASE)  # raiz del proyecto
    RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_PROYECTO, "datos", "crudo")
    RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")

    # Validar argumentos
    if len(sys.argv) < 2:
        print("[ERROR] Uso: python analisis_estructura_gen.py <genome_basename>")
        sys.exit(1)

    GENOME_BASENAME = sys.argv[1]
    ARCHIVO_GENBANK = os.path.join(RUTA_DATOS_CRUDO, f"{GENOME_BASENAME}.gb")
    
    # Crear directorio de resultados
    os.makedirs(RUTA_RESULTADOS, exist_ok=True)
    
    # Validar archivo antes de continuar
    if not os.path.exists(ARCHIVO_GENBANK):
        print(f"[ERROR] Archivo no encontrado: {ARCHIVO_GENBANK}")
        if os.path.exists(RUTA_DATOS_CRUDO):
            print(f"[INFO] Archivos disponibles: {os.listdir(RUTA_DATOS_CRUDO)}")
        sys.exit(1)

except Exception as e:
    print(f"[ERROR] Fallo durante inicializacion: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)


# =============================================================================
# FUNCIONES DE ANALISIS
# =============================================================================

def analizar_composicion_genomica(registro):
    """
    Calcula la composicion del genoma: bases CDS, tRNA, rRNA, etc.

    Returns:
        dict: Composicion con bases y porcentajes por tipo
    """
    print("\n" + "=" * 70)
    print("ANALISIS DE COMPOSICION GENOMICA")
    print("=" * 70)

    longitud_genoma = len(registro.seq)
    print(f"\n  Longitud del genoma: {longitud_genoma:,} pb")

    # Rastrear bases cubiertas por cada tipo de feature
    tipos = {
        "CDS": set(),
        "tRNA": set(),
        "rRNA": set(),
        "repeat_region": set(),
        "ncRNA": set(),
        "misc_RNA": set(),
        "tmRNA": set(),
    }

    for feature in registro.features:
        ft = feature.type
        if ft in tipos:
            start = int(feature.location.start)
            end = int(feature.location.end)
            tipos[ft].update(range(start, end))

    # Calcular bases unicas por tipo (evitar doble conteo)
    bases_cds = len(tipos["CDS"])
    bases_trna = len(tipos["tRNA"] - tipos["CDS"])
    bases_rrna = len(tipos["rRNA"] - tipos["CDS"] - tipos["tRNA"])
    bases_repeat = len(tipos["repeat_region"] - tipos["CDS"] - tipos["tRNA"] - tipos["rRNA"])
    bases_ncrna = len((tipos["ncRNA"] | tipos["misc_RNA"] | tipos["tmRNA"]) - tipos["CDS"] - tipos["tRNA"] - tipos["rRNA"])

    todas_anotadas = tipos["CDS"] | tipos["tRNA"] | tipos["rRNA"] | tipos["repeat_region"] | tipos["ncRNA"] | tipos["misc_RNA"] | tipos["tmRNA"]
    bases_intergenico = longitud_genoma - len(todas_anotadas)

    composicion = {
        "codificante_cds": {
            "bases": bases_cds,
            "porcentaje": round(bases_cds / longitud_genoma * 100, 2)
        },
        "tRNA": {
            "bases": bases_trna,
            "porcentaje": round(bases_trna / longitud_genoma * 100, 2)
        },
        "rRNA": {
            "bases": bases_rrna,
            "porcentaje": round(bases_rrna / longitud_genoma * 100, 2)
        },
        "regiones_repetitivas": {
            "bases": bases_repeat,
            "porcentaje": round(bases_repeat / longitud_genoma * 100, 2)
        },
        "otro_RNA_funcional": {
            "bases": bases_ncrna,
            "porcentaje": round(bases_ncrna / longitud_genoma * 100, 2)
        },
        "intergenico": {
            "bases": bases_intergenico,
            "porcentaje": round(bases_intergenico / longitud_genoma * 100, 2)
        }
    }

    print(f"\n  Composicion del genoma:")
    print(f"  {'Tipo':<30} {'Bases':>12} {'Porcentaje':>12}")
    print("  " + "-" * 54)
    for tipo, datos in composicion.items():
        print(f"  {tipo:<30} {datos['bases']:>12,} {datos['porcentaje']:>11.2f}%")

    return composicion


def contar_features(registro):
    """
    Cuenta todos los tipos de features en el GenBank.

    Returns:
        dict: Conteo por tipo de feature
    """
    conteo = Counter()
    for feature in registro.features:
        conteo[feature.type] += 1

    print(f"\n  Features anotados en el GenBank:")
    for tipo, n in conteo.most_common():
        if tipo != "source":
            print(f"    {tipo:<25} {n:>6}")

    gc.collect()  # Liberar memoria
    return dict(conteo)


def identificar_operones(registro):
    """
    Identifica operones putativos: grupos de CDS adyacentes en la misma hebra
    con distancias intergenicas menores a 50 pb.

    Returns:
        dict: Operones encontrados
    """
    print("\n" + "=" * 70)
    print("IDENTIFICACION DE OPERONES PUTATIVOS")
    print("=" * 70)

    # Extraer CDS ordenados por posicion
    genes_cds = []
    for feature in registro.features:
        if feature.type == "CDS":
            genes_cds.append({
                "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
                "gene": feature.qualifiers.get("gene", [""])[0],
                "product": feature.qualifiers.get("product", [""])[0],
                "start": int(feature.location.start),
                "end": int(feature.location.end),
                "strand": "+" if feature.location.strand == 1 else "-"
            })

    genes_cds.sort(key=lambda g: g["start"])

    # Buscar grupos de genes adyacentes en misma hebra con gap < 50bp
    UMBRAL_DISTANCIA = 50  # pb
    operones = []
    operon_actual = [genes_cds[0]] if genes_cds else []

    for i in range(1, len(genes_cds)):
        gen_prev = genes_cds[i - 1]
        gen_actual = genes_cds[i]

        misma_hebra = gen_prev["strand"] == gen_actual["strand"]
        distancia = gen_actual["start"] - gen_prev["end"]

        if misma_hebra and 0 <= distancia <= UMBRAL_DISTANCIA:
            operon_actual.append(gen_actual)
        else:
            if len(operon_actual) >= 2:
                operones.append(operon_actual)
            operon_actual = [gen_actual]

    # Ultimo operon
    if len(operon_actual) >= 2:
        operones.append(operon_actual)

    # Formatear resultados - solo guardar lo esencial para ahorrar memoria
    operones_resultado = []
    for i, operon in enumerate(operones):
        nombres_genes = [g["gene"] or g["locus_tag"] for g in operon]
        operones_resultado.append({
            "numero": i + 1,
            "genes": nombres_genes,
            "num_genes": len(operon),
            "hebra": operon[0]["strand"],
            "inicio": operon[0]["start"],
            "fin": operon[-1]["end"],
            "longitud_total": operon[-1]["end"] - operon[0]["start"]
        })
        
        if (i + 1) % 100 == 0:
            gc.collect()

    # Estadisticas
    total_operones = len(operones_resultado)
    genes_en_operones = sum(op["num_genes"] for op in operones_resultado)
    max_genes = max((op["num_genes"] for op in operones_resultado), default=0)

    print(f"\n  Operones putativos encontrados: {total_operones}")
    print(f"  Genes organizados en operones: {genes_en_operones} ({round(genes_en_operones / len(genes_cds) * 100, 1) if genes_cds else 0}%)")
    print(f"  Operon mas grande: {max_genes} genes")
    print(f"  Umbral de distancia: {UMBRAL_DISTANCIA} pb")

    # Mostrar top 10 operones mas grandes
    top_operones = sorted(operones_resultado, key=lambda x: x["num_genes"], reverse=True)[:10]
    print(f"\n  Top 10 operones mas grandes:")
    print(f"  {'#':<4} {'Genes':>6} {'Hebra':>6} {'Inicio':>12} {'Fin':>12} {'Nombres'}")
    print("  " + "-" * 70)
    for op in top_operones:
        nombres = ", ".join(op["genes"][:5])
        if len(op["genes"]) > 5:
            nombres += f" (+{len(op['genes']) - 5})"
        print(f"  {op['numero']:<4} {op['num_genes']:>6} {op['hebra']:>6} {op['inicio']:>12,} {op['fin']:>12,} {nombres}")

    gc.collect()  # Liberar memoria después de procesar operones
    return {
        "total": total_operones,
        "genes_en_operones": genes_en_operones,
        "porcentaje_genes_en_operones": round(genes_en_operones / len(genes_cds) * 100, 1) if genes_cds else 0,
        "operon_mas_grande": max_genes,
        "umbral_distancia_pb": UMBRAL_DISTANCIA,
        "operones": operones_resultado[:30],  # Limitar a 30 para ahorrar memoria
        "total_cds_analizados": len(genes_cds)
    }


def generar_datos_educativos():
    """
    Genera contenido educativo sobre estructura del gen bacteriano.

    Returns:
        dict: Datos educativos
    """
    return {
        "titulo": "Estructura del Gen Bacteriano vs Eucariota",
        "por_que_no_intrones": {
            "titulo": "Por que las bacterias no tienen intrones?",
            "explicacion": [
                "Las bacterias carecen del complejo de splicing (spliceosoma) necesario para procesar intrones. Este complejo de ribonucleoproteinas es exclusivo de celulas eucariotas.",
                "En bacterias, la transcripcion y traduccion ocurren simultaneamente (acoplamiento transcripcion-traduccion). No hay tiempo para eliminar intrones entre ambos procesos.",
                "La presion evolutiva favorece genomas compactos en procariotas. Los intrones representarian ADN 'extra' que consume energia sin beneficio directo, lo cual es desfavorable en organismos con alta tasa de replicacion.",
                "Los genes bacterianos se organizan frecuentemente en operones: multiples genes bajo un mismo promotor transcritos como un solo ARNm policistronico. Esta organizacion maximiza la eficiencia y es incompatible con el splicing de intrones.",
                "La densidad genica en bacterias es muy alta (85-90% del genoma codifica), a diferencia de eucariotas donde solo el 1-2% del genoma es codificante."
            ],
            "nota": "Excepcion: Algunos intrones del Grupo I y II se encuentran raramente en bacterias, especialmente en genes de ARNr y ARNt, pero son autocataliticos (se eliminan solos) y no requieren spliceosoma."
        },
        "estructura_gen_bacteriano": {
            "titulo": "Estructura de un Gen Bacteriano Tipico",
            "elementos": [
                {
                    "nombre": "Promotor",
                    "posicion": "Upstream (-35 y -10)",
                    "descripcion": "Region de reconocimiento por la ARN polimerasa. Contiene la caja -35 (TTGACA) y la caja -10 o Pribnow (TATAAT).",
                    "color": "#f59e0b"
                },
                {
                    "nombre": "RBS (Shine-Dalgarno)",
                    "posicion": "~8 pb antes del ATG",
                    "descripcion": "Sitio de union al ribosoma. Secuencia consenso AGGAGG complementaria al ARNr 16S. Esencial para el inicio de la traduccion.",
                    "color": "#8b5cf6"
                },
                {
                    "nombre": "Codon de Inicio (ATG)",
                    "posicion": "Inicio del CDS",
                    "descripcion": "Marca el inicio de la secuencia codificante. En bacterias, casi siempre es ATG (codifica formil-metionina, fMet).",
                    "color": "#10b981"
                },
                {
                    "nombre": "CDS (Secuencia Codificante)",
                    "posicion": "Cuerpo del gen",
                    "descripcion": "Region continua sin intrones que codifica la proteina. Cada 3 nucleotidos forman un codon que especifica un aminoacido.",
                    "color": "#3b82f6"
                },
                {
                    "nombre": "Codon de Parada",
                    "posicion": "Final del CDS",
                    "descripcion": "Senala el fin de la traduccion. En bacterias los mas comunes son TAA (mas frecuente), TGA y TAG.",
                    "color": "#ef4444"
                },
                {
                    "nombre": "Terminador",
                    "posicion": "Downstream del CDS",
                    "descripcion": "Estructura de horquilla (stem-loop) en el ARNm que causa la disociacion de la ARN polimerasa. Puede ser Rho-dependiente o Rho-independiente.",
                    "color": "#6b7280"
                }
            ]
        },
        "comparacion_eucariota": {
            "titulo": "Diferencias con Genes Eucariotas",
            "diferencias": [
                {"aspecto": "Intrones", "bacteria": "Ausentes (gen continuo)", "eucariota": "Presentes (requieren splicing)"},
                {"aspecto": "ARNm", "bacteria": "Policistronico (multiples genes)", "eucariota": "Monocistronico (un gen)"},
                {"aspecto": "Procesamiento", "bacteria": "Minimo (sin cap ni poly-A)", "eucariota": "Complejo (cap 5', poly-A 3', splicing)"},
                {"aspecto": "Transcripcion-Traduccion", "bacteria": "Simultaneas (acopladas)", "eucariota": "Separadas (nucleo vs citoplasma)"},
                {"aspecto": "Densidad genica", "bacteria": "Alta (~87% codificante)", "eucariota": "Baja (~1.5% codificante en humanos)"},
                {"aspecto": "Operones", "bacteria": "Frecuentes", "eucariota": "Raros (solo en nematodos)"},
                {"aspecto": "Tamano genoma", "bacteria": "Compacto (1-10 Mb)", "eucariota": "Grande (10 Mb - 150 Gb)"}
            ]
        }
    }


def analizar_gc_por_region(registro, composicion):
    """
    Calcula el contenido GC en regiones codificantes vs no codificantes.

    Returns:
        dict: GC por region
    """
    seq = str(registro.seq).upper()
    longitud = len(seq)

    # GC total
    gc_total = gc_fraction(registro.seq) * 100

    # GC en CDS
    bases_cds = []
    for feature in registro.features:
        if feature.type == "CDS":
            start = int(feature.location.start)
            end = int(feature.location.end)
            bases_cds.append(seq[start:end])

    seq_cds = "".join(bases_cds)
    gc_cds = (seq_cds.count("G") + seq_cds.count("C")) / len(seq_cds) * 100 if seq_cds else 0

    # GC en intergenico (aproximado)
    gc_intergenico = 0
    if composicion["intergenico"]["bases"] > 0:
        # Calcular GC excluyendo CDS
        total_gc = seq.count("G") + seq.count("C")
        gc_en_cds = seq_cds.count("G") + seq_cds.count("C")
        gc_intergenico_abs = total_gc - gc_en_cds
        bases_no_cds = longitud - len(seq_cds)
        gc_intergenico = gc_intergenico_abs / bases_no_cds * 100 if bases_no_cds > 0 else 0

    return {
        "gc_total": round(gc_total, 2),
        "gc_codificante": round(gc_cds, 2),
        "gc_no_codificante": round(gc_intergenico, 2),
        "diferencia_gc": round(gc_cds - gc_intergenico, 2)
    }


# =============================================================================
# FUNCION PRINCIPAL
# =============================================================================

# =============================================================================
# FUNCION PRINCIPAL
# =============================================================================

def main():
    """Ejecuta el analisis de estructura genomica."""
    try:
        print("\n" + "=" * 70)
        print(f"ANALISIS DE ESTRUCTURA GENOMICA")
        print(f"Genoma: {GENOME_BASENAME}")
        print("=" * 70)
        print(f"\n  Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"  Archivo: {ARCHIVO_GENBANK}")

        if not os.path.exists(ARCHIVO_GENBANK):
            print(f"\n[ERROR] No se encuentra el archivo: {ARCHIVO_GENBANK}")
            print(f"[INFO] Directorio de datos: {RUTA_DATOS_CRUDO}")
            print(f"[INFO] Archivos disponibles en directorio: {os.listdir(RUTA_DATOS_CRUDO) if os.path.exists(RUTA_DATOS_CRUDO) else 'Directorio no existe'}")
            sys.exit(1)

        # Leer GenBank
        print("\n[PASO 1/4] Leyendo archivo GenBank...")
        try:
            registro = SeqIO.read(ARCHIVO_GENBANK, "genbank")
        except Exception as e:
            print(f"[ERROR] Fallo al leer GenBank: {e}")
            raise
        
        organismo = registro.annotations.get("organism", GENOME_BASENAME.replace("_", " ").title())
        print(f"  Organismo: {organismo}")
        print(f"  Longitud: {len(registro.seq):,} pb")

        # Contar features
        print("\n[PASO 2/4] Contando features anotados...")
        try:
            features_conteo = contar_features(registro)
        except Exception as e:
            print(f"[ERROR] Fallo al contar features: {e}")
            raise

        # Analizar composicion
        print("\n[PASO 3/4] Analizando composicion genomica...")
        try:
            composicion = analizar_composicion_genomica(registro)
            gc.collect()  # Liberar memoria
        except Exception as e:
            print(f"[ERROR] Fallo en analisis de composicion: {e}")
            raise

        # GC por region
        try:
            gc_regiones = analizar_gc_por_region(registro, composicion)
            gc.collect()  # Liberar memoria
        except Exception as e:
            print(f"[ERROR] Fallo en analisis de GC: {e}")
            raise
            
        print(f"\n  Contenido GC por region:")
        print(f"    Total:          {gc_regiones['gc_total']:.2f}%")
        print(f"    Codificante:    {gc_regiones['gc_codificante']:.2f}%")
        print(f"    No codificante: {gc_regiones['gc_no_codificante']:.2f}%")
        print(f"    Diferencia:     {gc_regiones['diferencia_gc']:.2f}%")

        # Identificar operones (OPCIONAL - continua si falla)
        print("\n[PASO 4/4] Identificando operones putativos...")
        operones = None
        try:
            print("  [INICIO] Comenzando identificacion de operones...")
            operones = identificar_operones(registro)
            print("  [EXITO] Identificacion de operones completada")
            gc.collect()  # Liberar memoria
        except MemoryError as e:
            print(f"[WARN] Error de memoria en operones (continuando sin operones): {type(e).__name__}")
            # Crear estructura minima de operones si falla
            operones = {
                "total": 0,
                "genes_en_operones": 0,
                "porcentaje_genes_en_operones": 0,
                "operon_mas_grande": 0,
                "umbral_distancia_pb": 50,
                "operones": [],
                "total_cds_analizados": features_conteo.get("CDS", 0),
                "omitido_por_error": True,
                "error_tipo": "MemoryError"
            }
            gc.collect()  # Asegurar liberacion de memoria despues del error
        except Exception as e:
            print(f"[WARN] Fallo en identificacion de operones (continuando sin operones): {type(e).__name__}: {str(e)[:100]}")
            # Crear estructura minima de operones si falla
            operones = {
                "total": 0,
                "genes_en_operones": 0,
                "porcentaje_genes_en_operones": 0,
                "operon_mas_grande": 0,
                "umbral_distancia_pb": 50,
                "operones": [],
                "total_cds_analizados": features_conteo.get("CDS", 0),
                "omitido_por_error": True,
                "error_tipo": type(e).__name__
            }
            gc.collect()  # Asegurar liberacion de memoria despues del error
        
        # *** Liberar memoria del genoma entero despues de procesar ***
        print("\n[MEMORIA] Limpiando registro GenBank de memoria...")
        # Guardar datos que necesitamos luego (longitud) antes de eliminar el registro
        longitud_genoma = len(registro.seq)
        del registro
        gc.collect()

        # Datos educativos
        try:
            educativo = generar_datos_educativos()
        except Exception as e:
            print(f"[ERROR] Fallo al generar datos educativos: {e}")
            raise

        # Compilar resultados
        resultados = {
            "fecha_analisis": datetime.now().isoformat(),
            "organismo": organismo,
            "genoma_basename": GENOME_BASENAME,
            "longitud_genoma": longitud_genoma,
            "composicion_genomica": composicion,
            "features_anotados": features_conteo,
            "gc_por_region": gc_regiones,
            "operones_putativos": operones,
            "educativo": educativo,
            "estadisticas": {
                "total_cds": features_conteo.get("CDS", 0),
                "total_trna": features_conteo.get("tRNA", 0),
                "total_rrna": features_conteo.get("rRNA", 0),
                "total_repeat_region": features_conteo.get("repeat_region", 0),
                "total_ncrna": features_conteo.get("ncRNA", 0) + features_conteo.get("misc_RNA", 0),
                "densidad_codificante": composicion["codificante_cds"]["porcentaje"],
                "total_operones": operones["total"],
                "genes_en_operones": operones["genes_en_operones"]
            }
        }

        # Exportar con error handling
        print(f"\n[PASO 5/5] Exportando resultados a JSON...")
        try:
            os.makedirs(RUTA_RESULTADOS, exist_ok=True)
            archivo_salida = os.path.join(RUTA_RESULTADOS, f"analisis_estructura_{GENOME_BASENAME}.json")
            
            with open(archivo_salida, "w", encoding="utf-8") as f:
                json.dump(resultados, f, indent=2, ensure_ascii=False)
            
            print(f"  [OK] Archivo JSON creado: {archivo_salida}")
            file_size = os.path.getsize(archivo_salida)
            print(f"  [OK] Tamano del archivo: {file_size / 1024:.1f} KB")
        except Exception as e:
            print(f"[ERROR] Fallo al guardar JSON: {e}")
            raise

        print(f"\n" + "=" * 70)
        print(f"[OK] Analisis de estructura completado")
        print(f"     Archivo: analisis_estructura_{GENOME_BASENAME}.json")
        print("=" * 70 + "\n")
        
    except Exception as e:
        print(f"\n[ERROR] Excepcion durante analisis de estructura:")
        print(f"  {type(e).__name__}: {str(e)}")
        import traceback
        print("\n[TRACEBACK]")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

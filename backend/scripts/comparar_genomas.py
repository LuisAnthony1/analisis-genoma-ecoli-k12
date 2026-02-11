#!/usr/bin/env python3
"""
Comparacion de Genomas Bacterianos

Este script compara dos genomas bacterianos cualesquiera a partir de sus
archivos GenBank descargados.

Uso: python comparar_genomas.py <basename_1> <basename_2>
Ejemplo: python comparar_genomas.py escherichia_coli_str_k_12_substr_mg1655 salmonella_enterica_subsp_enterica_serovar_typhimu

Funcionalidades:
- Comparar metricas generales (tamano, GC, densidad genica)
- Analizar diferencias en numero y distribucion de genes
- Identificar caracteristicas de patogenicidad
- Calcular distancias intergenicas (indicador de islas genomicas)
- Exportar comparaciones a CSV y JSON

Fecha: 2026
Proyecto: Analisis comparativo de genomas bacterianos
"""

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os
import sys
import json
import csv
from datetime import datetime
from collections import Counter
import statistics


# CONFIGURACION
# =============================================================================

DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # backend/
DIRECTORIO_PROYECTO = os.path.dirname(DIRECTORIO_BASE)  # raiz del proyecto
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_PROYECTO, "datos", "crudo")
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")

# Parametros de linea de comandos
if len(sys.argv) < 3:
    print("[ERROR] Uso: python comparar_genomas.py <basename_1> <basename_2>")
    print("  Ejemplo: python comparar_genomas.py ecoli_k12 salmonella_lt2")
    sys.exit(1)

BASENAME_1 = sys.argv[1]
BASENAME_2 = sys.argv[2]

ARCHIVO_GB_1 = os.path.join(RUTA_DATOS_CRUDO, f"{BASENAME_1}.gb")
ARCHIVO_GB_2 = os.path.join(RUTA_DATOS_CRUDO, f"{BASENAME_2}.gb")

# Palabras clave para detectar genes de virulencia/patogenicidad
PALABRAS_VIRULENCIA = [
    "invasion", "invasin", "virulence", "pathogen", "toxin", "secretion",
    "type III", "T3SS", "SPI", "effector", "adhesin", "fimbri", "pili",
    "flagell", "motility", "iron", "siderophore", "enterobactin",
    "resistance", "antibiotic", "drug"
]


# FUNCIONES DE CARGA
# =============================================================================

def crear_directorios():
    """Crea los directorios de resultados si no existen."""
    if not os.path.exists(RUTA_RESULTADOS):
        os.makedirs(RUTA_RESULTADOS)
        print(f"[INFO] Directorio creado: {RUTA_RESULTADOS}")


def verificar_archivos():
    """
    Verifica que existan los archivos GenBank de ambos genomas.

    Returns:
        tuple: (bool, list) - Si ambos existen y lista de faltantes
    """
    faltantes = []

    if not os.path.exists(ARCHIVO_GB_1):
        faltantes.append(f"{BASENAME_1}.gb")
    if not os.path.exists(ARCHIVO_GB_2):
        faltantes.append(f"{BASENAME_2}.gb")

    return len(faltantes) == 0, faltantes


def cargar_genbank(archivo_gb, nombre):
    """
    Carga un archivo GenBank.

    Args:
        archivo_gb: Ruta al archivo GenBank
        nombre: Nombre para mostrar en logs

    Returns:
        SeqRecord: Registro del genoma
    """
    print(f"[INFO] Cargando {nombre}...")
    registro = SeqIO.read(archivo_gb, "genbank")
    print(f"       Longitud: {len(registro.seq):,} pares de bases")
    return registro


def obtener_nombre_organismo(registro, basename):
    """Obtiene el nombre del organismo desde el registro GenBank."""
    nombre = registro.annotations.get("organism", "")
    if nombre:
        return nombre
    return basename.replace("_", " ").title()


# FUNCIONES DE EXTRACCION
# =============================================================================

def extraer_genes_con_info(registro, nombre_organismo):
    """
    Extrae genes con informacion detallada para comparacion.

    Args:
        registro: SeqRecord del genoma
        nombre_organismo: Nombre del organismo

    Returns:
        list: Lista de genes con informacion
    """
    genes = []
    secuencia_genoma = registro.seq

    for feature in registro.features:
        if feature.type == "CDS":
            inicio = int(feature.location.start)
            fin = int(feature.location.end)
            longitud = fin - inicio

            secuencia_gen = feature.location.extract(secuencia_genoma)
            gc_gen = gc_fraction(secuencia_gen) * 100

            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
            nombre_gen = feature.qualifiers.get("gene", [""])[0]
            producto = feature.qualifiers.get("product", [""])[0]

            # Detectar si es gen de virulencia
            es_virulencia = any(
                palabra.lower() in producto.lower()
                for palabra in PALABRAS_VIRULENCIA
            )

            genes.append({
                "locus_tag": locus_tag,
                "nombre_gen": nombre_gen,
                "inicio": inicio,
                "fin": fin,
                "longitud_pb": longitud,
                "contenido_gc": round(gc_gen, 2),
                "producto": producto,
                "es_virulencia": es_virulencia,
                "organismo": nombre_organismo
            })

    return genes


def calcular_distancias_intergenicas(genes):
    """
    Calcula las distancias entre genes consecutivos.
    Distancias grandes pueden indicar islas genomicas o regiones adquiridas.

    Args:
        genes: Lista de genes ordenados por posicion

    Returns:
        dict: Estadisticas de distancias intergenicas
    """
    # Ordenar genes por posicion de inicio
    genes_ordenados = sorted(genes, key=lambda g: g["inicio"])

    distancias = []
    regiones_grandes = []  # Distancias > 5000 pb

    for i in range(1, len(genes_ordenados)):
        gen_anterior = genes_ordenados[i - 1]
        gen_actual = genes_ordenados[i]

        distancia = gen_actual["inicio"] - gen_anterior["fin"]

        if distancia > 0:  # Solo distancias positivas (genes no solapados)
            distancias.append(distancia)

            if distancia > 5000:  # Regiones intergenicas grandes
                regiones_grandes.append({
                    "gen_anterior": gen_anterior["locus_tag"],
                    "gen_siguiente": gen_actual["locus_tag"],
                    "distancia_pb": distancia,
                    "posicion_inicio": gen_anterior["fin"],
                    "posicion_fin": gen_actual["inicio"]
                })

    if not distancias:
        return {"error": "No se pudieron calcular distancias"}

    return {
        "promedio_pb": round(statistics.mean(distancias), 2),
        "mediana_pb": round(statistics.median(distancias), 2),
        "minima_pb": min(distancias),
        "maxima_pb": max(distancias),
        "desviacion_std": round(statistics.stdev(distancias), 2) if len(distancias) > 1 else 0,
        "total_regiones_intergenicas": len(distancias),
        "regiones_grandes_5kb": len(regiones_grandes),
        "detalle_regiones_grandes": regiones_grandes[:10]  # Top 10
    }


def contar_genes_virulencia(genes):
    """
    Cuenta y lista genes relacionados con virulencia.

    Args:
        genes: Lista de genes

    Returns:
        dict: Conteo y lista de genes de virulencia
    """
    genes_virulencia = [g for g in genes if g["es_virulencia"]]

    # Agrupar por categoria
    categorias = Counter()
    for gen in genes_virulencia:
        producto_lower = gen["producto"].lower()
        if "secretion" in producto_lower or "t3ss" in producto_lower:
            categorias["Sistema de secrecion"] += 1
        elif "invasion" in producto_lower or "invasin" in producto_lower:
            categorias["Invasion celular"] += 1
        elif "toxin" in producto_lower:
            categorias["Toxinas"] += 1
        elif "fimbri" in producto_lower or "pili" in producto_lower or "adhesin" in producto_lower:
            categorias["Adhesion/Fimbrias"] += 1
        elif "flagell" in producto_lower or "motility" in producto_lower:
            categorias["Motilidad/Flagelos"] += 1
        elif "iron" in producto_lower or "siderophore" in producto_lower:
            categorias["Captacion de hierro"] += 1
        elif "resistance" in producto_lower or "antibiotic" in producto_lower:
            categorias["Resistencia"] += 1
        else:
            categorias["Otros factores de virulencia"] += 1

    return {
        "total": len(genes_virulencia),
        "categorias": dict(categorias),
        "ejemplos": [
            {"locus": g["locus_tag"], "gen": g["nombre_gen"], "producto": g["producto"]}
            for g in genes_virulencia[:15]
        ]
    }


# FUNCIONES DE COMPARACION
# =============================================================================

def comparar_metricas_generales(registro_1, registro_2, genes_1, genes_2, nombre_1, nombre_2):
    """
    Compara metricas generales entre ambos genomas.

    Returns:
        dict: Comparacion de metricas
    """
    print("\n" + "=" * 70)
    print("COMPARACION DE METRICAS GENERALES")
    print("=" * 70)

    # Calcular metricas para genoma 1
    long_1 = len(registro_1.seq)
    gc_1 = gc_fraction(registro_1.seq) * 100
    total_cds_1 = len(genes_1)
    pb_codificante_1 = sum(g["longitud_pb"] for g in genes_1)
    densidad_1 = (pb_codificante_1 / long_1) * 100 if long_1 > 0 else 0

    # Calcular metricas para genoma 2
    long_2 = len(registro_2.seq)
    gc_2 = gc_fraction(registro_2.seq) * 100
    total_cds_2 = len(genes_2)
    pb_codificante_2 = sum(g["longitud_pb"] for g in genes_2)
    densidad_2 = (pb_codificante_2 / long_2) * 100 if long_2 > 0 else 0

    # Calcular diferencias
    dif_longitud = long_2 - long_1
    dif_gc = gc_2 - gc_1
    dif_genes = total_cds_2 - total_cds_1
    dif_densidad = densidad_2 - densidad_1

    # Usar "ecoli" y "salmonella" como claves para compatibilidad con dashboard
    comparacion = {
        "ecoli": {
            "nombre": nombre_1,
            "longitud_genoma_pb": long_1,
            "contenido_gc_porcentaje": round(gc_1, 2),
            "total_genes_cds": total_cds_1,
            "pb_codificante": pb_codificante_1,
            "densidad_genica_porcentaje": round(densidad_1, 2),
            "tamano_promedio_gen_pb": round(statistics.mean([g["longitud_pb"] for g in genes_1]), 2) if genes_1 else 0
        },
        "salmonella": {
            "nombre": nombre_2,
            "longitud_genoma_pb": long_2,
            "contenido_gc_porcentaje": round(gc_2, 2),
            "total_genes_cds": total_cds_2,
            "pb_codificante": pb_codificante_2,
            "densidad_genica_porcentaje": round(densidad_2, 2),
            "tamano_promedio_gen_pb": round(statistics.mean([g["longitud_pb"] for g in genes_2]), 2) if genes_2 else 0
        },
        "diferencias": {
            "longitud_pb": dif_longitud,
            "longitud_porcentaje": round((dif_longitud / long_1) * 100, 2) if long_1 > 0 else 0,
            "contenido_gc": round(dif_gc, 2),
            "total_genes": dif_genes,
            "densidad_genica": round(dif_densidad, 2)
        }
    }

    # Mostrar resultados
    n1_short = nombre_1[:20]
    n2_short = nombre_2[:20]
    print(f"\n  {'Metrica':<35} {n1_short:>20} {n2_short:>20} {'Diferencia':>15}")
    print("  " + "-" * 90)

    print(f"  {'Longitud del genoma (pb)':<35} {long_1:>20,} {long_2:>20,} {dif_longitud:>+15,}")
    print(f"  {'Contenido GC (%)':<35} {gc_1:>20.2f} {gc_2:>20.2f} {dif_gc:>+15.2f}")
    print(f"  {'Total genes (CDS)':<35} {total_cds_1:>20,} {total_cds_2:>20,} {dif_genes:>+15,}")
    print(f"  {'Densidad genica (%)':<35} {densidad_1:>20.2f} {densidad_2:>20.2f} {dif_densidad:>+15.2f}")

    print(f"\n  [INTERPRETACION]")
    print(f"  - Genoma 2 tiene {abs(dif_longitud):,} pb ({abs(comparacion['diferencias']['longitud_porcentaje']):.1f}%) "
          f"{'mas' if dif_longitud > 0 else 'menos'} que Genoma 1")
    print(f"  - Genoma 2 tiene {abs(dif_genes):,} genes {'mas' if dif_genes > 0 else 'menos'} que Genoma 1")
    print(f"  - El contenido GC difiere en {abs(dif_gc):.2f}%")

    return comparacion


def comparar_virulencia(genes_1, genes_2, nombre_1, nombre_2):
    """
    Compara genes de virulencia entre ambos organismos.

    Returns:
        dict: Comparacion de genes de virulencia
    """
    print("\n" + "=" * 70)
    print("COMPARACION DE GENES DE VIRULENCIA")
    print("=" * 70)
    print("  (Genes relacionados con patogenicidad, invasion, toxinas, etc.)")

    virulencia_1 = contar_genes_virulencia(genes_1)
    virulencia_2 = contar_genes_virulencia(genes_2)

    comparacion = {
        "ecoli": virulencia_1,
        "salmonella": virulencia_2,
        "diferencia_total": virulencia_2["total"] - virulencia_1["total"]
    }

    print(f"\n  CONTEO TOTAL DE GENES DE VIRULENCIA:")
    print(f"    {nombre_1}: {virulencia_1['total']:,} genes")
    print(f"    {nombre_2}: {virulencia_2['total']:,} genes")
    print(f"    Diferencia: {comparacion['diferencia_total']:+,} genes")

    print(f"\n  DESGLOSE POR CATEGORIA EN {nombre_2}:")
    for categoria, cantidad in sorted(virulencia_2["categorias"].items(), key=lambda x: -x[1]):
        cant_1 = virulencia_1["categorias"].get(categoria, 0)
        print(f"    {categoria:<35} {cantidad:>5} ({nombre_1[:15]}: {cant_1})")

    if virulencia_2["ejemplos"]:
        print(f"\n  EJEMPLOS DE GENES DE VIRULENCIA EN {nombre_2} (top 10):")
        for i, gen in enumerate(virulencia_2["ejemplos"][:10], 1):
            print(f"    {i:2}. {gen['locus']}: {gen['producto'][:60]}")

    return comparacion


def comparar_distancias_intergenicas(genes_1, genes_2, nombre_1, nombre_2):
    """
    Compara las distancias intergenicas entre ambos genomas.
    Distancias grandes pueden indicar islas de patogenicidad.

    Returns:
        dict: Comparacion de distancias intergenicas
    """
    print("\n" + "=" * 70)
    print("COMPARACION DE DISTANCIAS INTERGENICAS")
    print("=" * 70)
    print("  (Las distancias grandes pueden indicar islas genomicas/patogenicidad)")

    dist_1 = calcular_distancias_intergenicas(genes_1)
    dist_2 = calcular_distancias_intergenicas(genes_2)

    comparacion = {
        "ecoli": dist_1,
        "salmonella": dist_2
    }

    n1_short = nombre_1[:20]
    n2_short = nombre_2[:20]
    print(f"\n  {'Metrica':<40} {n1_short:>20} {n2_short:>20}")
    print("  " + "-" * 80)
    print(f"  {'Distancia intergenica promedio (pb)':<40} {dist_1['promedio_pb']:>20.2f} {dist_2['promedio_pb']:>20.2f}")
    print(f"  {'Distancia intergenica mediana (pb)':<40} {dist_1['mediana_pb']:>20.2f} {dist_2['mediana_pb']:>20.2f}")
    print(f"  {'Distancia maxima (pb)':<40} {dist_1['maxima_pb']:>20,} {dist_2['maxima_pb']:>20,}")
    print(f"  {'Regiones intergenicas > 5 kb':<40} {dist_1['regiones_grandes_5kb']:>20} {dist_2['regiones_grandes_5kb']:>20}")

    return comparacion


def comparar_uso_codones(registro_1, registro_2, genes_1, genes_2, nombre_1, nombre_2):
    """
    Compara el uso de codones entre ambos genomas.

    Returns:
        dict: Comparacion de uso de codones
    """
    print("\n" + "=" * 70)
    print("COMPARACION DE USO DE CODONES")
    print("=" * 70)

    def contar_codones(registro, genes):
        """Cuenta codones en todas las CDS del genoma."""
        conteo = Counter()
        secuencia = registro.seq
        for gen in genes:
            try:
                seq_gen = str(secuencia[gen["inicio"]:gen["fin"]])
                for i in range(0, len(seq_gen) - 2, 3):
                    codon = seq_gen[i:i+3].upper()
                    if len(codon) == 3 and all(c in "ATGC" for c in codon):
                        conteo[codon] += 1
            except Exception:
                continue
        return conteo

    codones_1 = contar_codones(registro_1, genes_1)
    codones_2 = contar_codones(registro_2, genes_2)

    total_1 = sum(codones_1.values())
    total_2 = sum(codones_2.values())

    # Top 10 mas usados en cada genoma
    top10_1 = codones_1.most_common(10)
    top10_2 = codones_2.most_common(10)

    # Frecuencias relativas
    freq_1 = {c: round(n / total_1 * 100, 3) if total_1 > 0 else 0 for c, n in codones_1.items()}
    freq_2 = {c: round(n / total_2 * 100, 3) if total_2 > 0 else 0 for c, n in codones_2.items()}

    # Codones con mayor diferencia
    todos_codones = set(list(freq_1.keys()) + list(freq_2.keys()))
    diferencias = []
    for codon in todos_codones:
        f1 = freq_1.get(codon, 0)
        f2 = freq_2.get(codon, 0)
        diferencias.append({
            "codon": codon,
            "frecuencia_1": f1,
            "frecuencia_2": f2,
            "diferencia_abs": round(abs(f2 - f1), 3)
        })
    diferencias.sort(key=lambda x: x["diferencia_abs"], reverse=True)

    print(f"\n  Top 10 codones mas usados:")
    print(f"  {'Codon':<10} {nombre_1[:15]:>15} {nombre_2[:15]:>15}")
    all_codons = set([c for c, _ in top10_1] + [c for c, _ in top10_2])
    for codon in list(all_codons)[:10]:
        f1 = freq_1.get(codon, 0)
        f2 = freq_2.get(codon, 0)
        print(f"  {codon:<10} {f1:>14.3f}% {f2:>14.3f}%")

    return {
        "ecoli": {
            "total_codones": total_1,
            "top_10": [{"codon": c, "conteo": n, "frecuencia": freq_1.get(c, 0)} for c, n in top10_1]
        },
        "salmonella": {
            "total_codones": total_2,
            "top_10": [{"codon": c, "conteo": n, "frecuencia": freq_2.get(c, 0)} for c, n in top10_2]
        },
        "mayores_diferencias": diferencias[:15]
    }


def comparar_distribucion_tamanos(genes_1, genes_2, nombre_1, nombre_2):
    """
    Compara la distribucion de tamanos de genes entre genomas.

    Returns:
        dict: Estadisticas de tamanos de genes
    """
    print("\n" + "=" * 70)
    print("COMPARACION DE DISTRIBUCION DE TAMANOS DE GENES")
    print("=" * 70)

    def calc_stats(genes):
        tamanos = [g["longitud_pb"] for g in genes]
        if not tamanos:
            return {}
        return {
            "promedio": round(statistics.mean(tamanos), 1),
            "mediana": round(statistics.median(tamanos), 1),
            "minimo": min(tamanos),
            "maximo": max(tamanos),
            "desviacion_std": round(statistics.stdev(tamanos), 1) if len(tamanos) > 1 else 0,
            "q25": round(sorted(tamanos)[len(tamanos) // 4], 1),
            "q75": round(sorted(tamanos)[3 * len(tamanos) // 4], 1),
        }

    def calc_rangos(genes):
        rangos = {"<300": 0, "300-600": 0, "600-900": 0, "900-1200": 0, "1200-1500": 0, ">1500": 0}
        for g in genes:
            t = g["longitud_pb"]
            if t < 300: rangos["<300"] += 1
            elif t < 600: rangos["300-600"] += 1
            elif t < 900: rangos["600-900"] += 1
            elif t < 1200: rangos["900-1200"] += 1
            elif t < 1500: rangos["1200-1500"] += 1
            else: rangos[">1500"] += 1
        return rangos

    stats_1 = calc_stats(genes_1)
    stats_2 = calc_stats(genes_2)
    rangos_1 = calc_rangos(genes_1)
    rangos_2 = calc_rangos(genes_2)

    print(f"\n  {'Estadistica':<25} {nombre_1[:15]:>15} {nombre_2[:15]:>15}")
    print("  " + "-" * 55)
    for key in ["promedio", "mediana", "minimo", "maximo"]:
        print(f"  {key.title():<25} {stats_1.get(key, 0):>15} {stats_2.get(key, 0):>15}")

    return {
        "ecoli": {"estadisticas": stats_1, "rangos": rangos_1},
        "salmonella": {"estadisticas": stats_2, "rangos": rangos_2}
    }


def comparar_distribucion_gc(genes_1, genes_2, nombre_1, nombre_2):
    """
    Compara la distribucion de contenido GC por gen.

    Returns:
        dict: Distribucion de GC por gen
    """
    print("\n" + "=" * 70)
    print("COMPARACION DE DISTRIBUCION GC POR GEN")
    print("=" * 70)

    def calc_gc_stats(genes):
        gcs = [g["contenido_gc"] for g in genes if g.get("contenido_gc")]
        if not gcs:
            return {}
        return {
            "promedio": round(statistics.mean(gcs), 2),
            "mediana": round(statistics.median(gcs), 2),
            "minimo": round(min(gcs), 2),
            "maximo": round(max(gcs), 2),
            "desviacion_std": round(statistics.stdev(gcs), 2) if len(gcs) > 1 else 0,
        }

    def calc_rangos_gc(genes):
        rangos = {"<40%": 0, "40-45%": 0, "45-50%": 0, "50-55%": 0, "55-60%": 0, ">60%": 0}
        for g in genes:
            gc = g.get("contenido_gc", 0)
            if gc < 40: rangos["<40%"] += 1
            elif gc < 45: rangos["40-45%"] += 1
            elif gc < 50: rangos["45-50%"] += 1
            elif gc < 55: rangos["50-55%"] += 1
            elif gc < 60: rangos["55-60%"] += 1
            else: rangos[">60%"] += 1
        return rangos

    stats_1 = calc_gc_stats(genes_1)
    stats_2 = calc_gc_stats(genes_2)
    rangos_1 = calc_rangos_gc(genes_1)
    rangos_2 = calc_rangos_gc(genes_2)

    print(f"\n  {'Estadistica':<25} {nombre_1[:15]:>15} {nombre_2[:15]:>15}")
    print("  " + "-" * 55)
    for key in ["promedio", "mediana", "minimo", "maximo"]:
        print(f"  {key.title():<25} {stats_1.get(key, 0):>14.2f}% {stats_2.get(key, 0):>14.2f}%")

    return {
        "ecoli": {"estadisticas": stats_1, "rangos": rangos_1},
        "salmonella": {"estadisticas": stats_2, "rangos": rangos_2}
    }


def obtener_interpretacion_ia(resultados_completos, nombre_1, nombre_2):
    """
    Obtiene una interpretacion biologica de los resultados via Gemini.

    Returns:
        str: Interpretacion en texto o None
    """
    print("\n[INFO] Consultando IA para interpretacion biologica...")
    try:
        sys.path.insert(0, DIRECTORIO_PROYECTO)
        from backend.gemini_client import consultar_gemini

        # Preparar resumen para la IA
        metricas = resultados_completos.get("metricas_generales", {})
        ec = metricas.get("ecoli", {})
        sal = metricas.get("salmonella", {})
        vir = resultados_completos.get("genes_virulencia", {})

        resumen = f"""Comparacion genomica entre {nombre_1} y {nombre_2}:
- {nombre_1}: {ec.get('longitud_genoma_pb', 0):,} pb, {ec.get('total_genes_cds', 0)} CDS, GC {ec.get('contenido_gc_porcentaje', 0)}%, densidad {ec.get('densidad_genica_porcentaje', 0)}%
- {nombre_2}: {sal.get('longitud_genoma_pb', 0):,} pb, {sal.get('total_genes_cds', 0)} CDS, GC {sal.get('contenido_gc_porcentaje', 0)}%, densidad {sal.get('densidad_genica_porcentaje', 0)}%
- Genes virulencia: {nombre_1}={vir.get('ecoli', {}).get('total', 0)}, {nombre_2}={vir.get('salmonella', {}).get('total', 0)}"""

        prompt = f"""Eres un experto en genomica comparativa. Analiza esta comparacion y da una interpretacion biologica concisa en espanol (3-5 parrafos).

{resumen}

Incluye:
1. Interpretacion de las diferencias en tamano y contenido genomico
2. Significado evolutivo de las diferencias en GC%
3. Implicaciones de los genes de virulencia
4. Conclusiones biologicas principales

Responde de forma directa sin introducciones."""

        respuesta, error = consultar_gemini(prompt)
        if respuesta:
            print("[INFO] Interpretacion IA obtenida exitosamente")
            return respuesta
        else:
            print(f"[INFO] No se pudo obtener interpretacion IA: {error}")
            return None
    except Exception as e:
        print(f"[INFO] Error al consultar IA: {e}")
        return None


def generar_resumen_comparativo(metricas, virulencia, distancias, nombre_1, nombre_2):
    """
    Genera un resumen narrativo de la comparacion.

    Returns:
        dict: Resumen con conclusiones principales
    """
    print("\n" + "=" * 70)
    print(f"RESUMEN COMPARATIVO: {nombre_1} vs {nombre_2}")
    print("=" * 70)

    resumen = {
        "titulo": f"Comparacion genomica: {nombre_1} vs {nombre_2}",
        "diferencias_clave": []
    }

    # Diferencia 1: Tamano del genoma
    dif_tamano = metricas["diferencias"]["longitud_pb"]
    resumen["diferencias_clave"].append({
        "aspecto": "Tamano del genoma",
        "observacion": f"Diferencia de {abs(dif_tamano):,} pb",
        "significado": "ADN adicional puede contener genes de adaptacion"
    })

    # Diferencia 2: Genes de virulencia
    dif_virulencia = virulencia["diferencia_total"]
    resumen["diferencias_clave"].append({
        "aspecto": "Genes de virulencia",
        "observacion": f"Diferencia de {abs(dif_virulencia)} genes de virulencia",
        "significado": "Genes para invasion celular, toxinas, evasion inmune"
    })

    # Diferencia 3: Regiones intergenicas
    dist_1 = distancias["ecoli"]
    dist_2 = distancias["salmonella"]
    resumen["diferencias_clave"].append({
        "aspecto": "Regiones intergenicas grandes",
        "observacion": f"Genoma 1: {dist_1.get('regiones_grandes_5kb', 0)}, Genoma 2: {dist_2.get('regiones_grandes_5kb', 0)} regiones >5kb",
        "significado": "Regiones grandes pueden indicar islas genomicas adquiridas"
    })

    print(f"\n  DIFERENCIAS CLAVE:")
    for i, dif in enumerate(resumen["diferencias_clave"], 1):
        print(f"\n    {i}. {dif['aspecto'].upper()}")
        print(f"       Observacion: {dif['observacion']}")
        print(f"       Significado: {dif['significado']}")

    resumen["conclusiones"] = [
        "La comparacion genomica revela diferencias en contenido genico y organizacion",
        "Las diferencias en genes de virulencia indican distintas capacidades patogenicas",
        "Las regiones intergenicas grandes pueden indicar adquisicion horizontal de genes"
    ]

    print(f"\n  CONCLUSIONES:")
    for i, conclusion in enumerate(resumen["conclusiones"], 1):
        print(f"    {i}. {conclusion}")

    return resumen


# FUNCIONES DE EXPORTACION
# =============================================================================

def exportar_comparacion_json(datos, nombre_archivo):
    """Exporta la comparacion completa a JSON."""
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.json")

    with open(ruta, 'w', encoding='utf-8') as f:
        json.dump(datos, f, indent=2, ensure_ascii=False)

    print(f"       [OK] Guardado: {nombre_archivo}.json")


def exportar_comparacion_csv(metricas, nombre_archivo, nombre_1, nombre_2):
    """Exporta las metricas comparativas a CSV."""
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")

    with open(ruta, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Metrica', nombre_1, nombre_2, 'Diferencia', 'Unidad'])

        ec = metricas["ecoli"]
        sal = metricas["salmonella"]
        dif = metricas["diferencias"]

        filas = [
            ('Longitud_genoma', ec['longitud_genoma_pb'], sal['longitud_genoma_pb'], dif['longitud_pb'], 'pb'),
            ('Contenido_GC', ec['contenido_gc_porcentaje'], sal['contenido_gc_porcentaje'], dif['contenido_gc'], '%'),
            ('Total_genes_CDS', ec['total_genes_cds'], sal['total_genes_cds'], dif['total_genes'], 'genes'),
            ('Densidad_genica', ec['densidad_genica_porcentaje'], sal['densidad_genica_porcentaje'], dif['densidad_genica'], '%'),
            ('Tamano_promedio_gen', ec['tamano_promedio_gen_pb'], sal['tamano_promedio_gen_pb'],
             sal['tamano_promedio_gen_pb'] - ec['tamano_promedio_gen_pb'], 'pb'),
        ]

        for fila in filas:
            writer.writerow(fila)

    print(f"       [OK] Guardado: {nombre_archivo}.csv")


def exportar_genes_virulencia_csv(virulencia_datos, nombre_archivo):
    """Exporta lista de genes de virulencia a CSV."""
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")

    with open(ruta, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Locus_tag', 'Nombre_gen', 'Producto'])

        for gen in virulencia_datos["ejemplos"]:
            writer.writerow([gen['locus'], gen['gen'], gen['producto']])

    print(f"       [OK] Guardado: {nombre_archivo}.csv")


# FUNCION PRINCIPAL
# =============================================================================

def main():
    """Funcion principal que ejecuta la comparacion de genomas."""

    print("\n" + "=" * 70)
    print("COMPARACION DE GENOMAS BACTERIANOS")
    print("=" * 70)
    print(f"\n  Fecha de analisis: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Genoma 1: {BASENAME_1}")
    print(f"  Genoma 2: {BASENAME_2}")
    print("")

    # Crear directorios
    crear_directorios()

    # Verificar archivos
    archivos_ok, faltantes = verificar_archivos()
    if not archivos_ok:
        print("[ERROR] Faltan archivos GenBank:")
        for f in faltantes:
            print(f"        - {f}")
        print("\n[INFO] Descarga primero los genomas desde la interfaz web")
        return

    # Cargar genomas
    print("\n[PASO 1/4] Cargando archivos GenBank...")
    registro_1 = cargar_genbank(ARCHIVO_GB_1, BASENAME_1)
    registro_2 = cargar_genbank(ARCHIVO_GB_2, BASENAME_2)

    nombre_1 = obtener_nombre_organismo(registro_1, BASENAME_1)
    nombre_2 = obtener_nombre_organismo(registro_2, BASENAME_2)
    print(f"       Organismo 1: {nombre_1}")
    print(f"       Organismo 2: {nombre_2}")

    # Extraer genes
    print("\n[PASO 2/4] Extrayendo genes de ambos genomas...")
    genes_1 = extraer_genes_con_info(registro_1, nombre_1)
    genes_2 = extraer_genes_con_info(registro_2, nombre_2)
    print(f"       {nombre_1}: {len(genes_1):,} genes CDS extraidos")
    print(f"       {nombre_2}: {len(genes_2):,} genes CDS extraidos")

    # Realizar comparaciones
    print("\n[PASO 3/4] Realizando comparaciones...")

    metricas = comparar_metricas_generales(
        registro_1, registro_2,
        genes_1, genes_2,
        nombre_1, nombre_2
    )

    virulencia = comparar_virulencia(genes_1, genes_2, nombre_1, nombre_2)

    distancias = comparar_distancias_intergenicas(genes_1, genes_2, nombre_1, nombre_2)

    # Nuevas comparaciones
    uso_codones = comparar_uso_codones(registro_1, registro_2, genes_1, genes_2, nombre_1, nombre_2)
    dist_tamanos = comparar_distribucion_tamanos(genes_1, genes_2, nombre_1, nombre_2)
    dist_gc = comparar_distribucion_gc(genes_1, genes_2, nombre_1, nombre_2)

    resumen = generar_resumen_comparativo(metricas, virulencia, distancias, nombre_1, nombre_2)

    # Compilar resultados
    resultados_completos = {
        "fecha_analisis": datetime.now().isoformat(),
        "organismos_comparados": {
            "organismo_1": {"nombre": nombre_1, "basename": BASENAME_1},
            "organismo_2": {"nombre": nombre_2, "basename": BASENAME_2}
        },
        "metricas_generales": metricas,
        "genes_virulencia": virulencia,
        "distancias_intergenicas": distancias,
        "uso_codones": uso_codones,
        "distribucion_tamanos": dist_tamanos,
        "distribucion_gc": dist_gc,
        "resumen_interpretativo": resumen
    }

    # Obtener interpretacion IA
    interpretacion = obtener_interpretacion_ia(resultados_completos, nombre_1, nombre_2)
    if interpretacion:
        resultados_completos["interpretacion_ia"] = interpretacion

    # Nombre del archivo de salida
    nombre_comparacion = f"comparacion_{BASENAME_1}_vs_{BASENAME_2}"

    # Exportar resultados
    print("\n[PASO 4/4] Exportando resultados...")
    print("  (Los archivos se guardan en resultados/tablas/)\n")

    print("  [1/3] Exportando comparacion completa (JSON)...")
    exportar_comparacion_json(resultados_completos, nombre_comparacion)

    print("  [2/3] Exportando metricas comparativas (CSV)...")
    exportar_comparacion_csv(metricas, f"metricas_{nombre_comparacion}", nombre_1, nombre_2)

    print("  [3/3] Exportando genes de virulencia (CSV)...")
    exportar_genes_virulencia_csv(virulencia["salmonella"], f"genes_virulencia_{BASENAME_2}")

    # Resumen final
    print("\n" + "=" * 70)
    print("ARCHIVOS GENERADOS")
    print("=" * 70)
    print(f"    - {nombre_comparacion}.json   (comparacion completa)")
    print(f"    - metricas_{nombre_comparacion}.csv   (metricas lado a lado)")
    print(f"    - genes_virulencia_{BASENAME_2}.csv   (genes de virulencia)")

    print("\n" + "=" * 70)
    print("[OK] Comparacion de genomas completada exitosamente")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()

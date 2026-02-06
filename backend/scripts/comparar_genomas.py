#!/usr/bin/env python3
"""
Comparacion de Genomas Bacterianos: E. coli vs Salmonella

Este script compara los genomas de dos bacterias relacionadas:
1. Escherichia coli K-12 MG1655 (no patogena, cepa de laboratorio)
2. Salmonella enterica serovar Typhimurium LT2 (patogena, causa gastroenteritis)

Objetivo:
Demostrar como pequenas diferencias geneticas pueden crear patogenos.
E. coli y Salmonella divergieron hace ~100-160 millones de anos, pero comparten
~70% de sus genes. Las diferencias incluyen islas de patogenicidad, genes de
virulencia y factores de colonizacion.

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
import json
import csv
from datetime import datetime
from collections import Counter
import statistics


# CONFIGURACION
# =============================================================================

DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_BASE, "datos", "crudo")
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")

# Configuracion de organismos
ORGANISMOS = {
    "ecoli": {
        "nombre": "Escherichia coli K-12 MG1655",
        "nombre_corto": "ecoli_k12",
        "archivo_genbank": "ecoli_k12.gb",
        "archivo_analisis": "analisis_genes_ecoli_k12.json",
        "patogeno": False,
        "descripcion": "Cepa de laboratorio no patogena",
        "habitat": "Intestino de mamiferos (comensal)",
    },
    "salmonella": {
        "nombre": "Salmonella enterica serovar Typhimurium LT2",
        "nombre_corto": "salmonella_lt2",
        "archivo_genbank": "salmonella_lt2.gb",
        "archivo_analisis": "analisis_genes_salmonella_lt2.json",
        "patogeno": True,
        "descripcion": "Patogeno que causa gastroenteritis",
        "habitat": "Intestino - invasivo (patogeno)",
    }
}

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
    Verifica que existan los archivos GenBank de ambos organismos.

    Returns:
        tuple: (bool, list) - Si ambos existen y lista de faltantes
    """
    faltantes = []

    for key, org in ORGANISMOS.items():
        ruta = os.path.join(RUTA_DATOS_CRUDO, org["archivo_genbank"])
        if not os.path.exists(ruta):
            faltantes.append(org["nombre"])

    return len(faltantes) == 0, faltantes


def cargar_genbank(nombre_organismo):
    """
    Carga un archivo GenBank.

    Args:
        nombre_organismo: "ecoli" o "salmonella"

    Returns:
        SeqRecord: Registro del genoma
    """
    config = ORGANISMOS[nombre_organismo]
    ruta = os.path.join(RUTA_DATOS_CRUDO, config["archivo_genbank"])

    print(f"[INFO] Cargando {config['nombre']}...")
    registro = SeqIO.read(ruta, "genbank")
    print(f"       Longitud: {len(registro.seq):,} pares de bases")

    return registro


def cargar_analisis_previo(nombre_organismo):
    """
    Carga el archivo JSON de analisis previo si existe.

    Args:
        nombre_organismo: "ecoli" o "salmonella"

    Returns:
        dict o None: Datos del analisis previo
    """
    config = ORGANISMOS[nombre_organismo]
    ruta = os.path.join(RUTA_RESULTADOS, config["archivo_analisis"])

    if os.path.exists(ruta):
        with open(ruta, 'r', encoding='utf-8') as f:
            return json.load(f)
    return None


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
        "desviacion_std": round(statistics.stdev(distancias), 2),
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

def comparar_metricas_generales(registro_ecoli, registro_salmonella, genes_ecoli, genes_salmonella):
    """
    Compara metricas generales entre ambos genomas.

    Returns:
        dict: Comparacion de metricas
    """
    print("\n" + "=" * 70)
    print("COMPARACION DE METRICAS GENERALES")
    print("=" * 70)

    # Calcular metricas para E. coli
    long_ecoli = len(registro_ecoli.seq)
    gc_ecoli = gc_fraction(registro_ecoli.seq) * 100
    total_cds_ecoli = len(genes_ecoli)
    pb_codificante_ecoli = sum(g["longitud_pb"] for g in genes_ecoli)
    densidad_ecoli = (pb_codificante_ecoli / long_ecoli) * 100

    # Calcular metricas para Salmonella
    long_salmonella = len(registro_salmonella.seq)
    gc_salmonella = gc_fraction(registro_salmonella.seq) * 100
    total_cds_salmonella = len(genes_salmonella)
    pb_codificante_salmonella = sum(g["longitud_pb"] for g in genes_salmonella)
    densidad_salmonella = (pb_codificante_salmonella / long_salmonella) * 100

    # Calcular diferencias
    dif_longitud = long_salmonella - long_ecoli
    dif_gc = gc_salmonella - gc_ecoli
    dif_genes = total_cds_salmonella - total_cds_ecoli
    dif_densidad = densidad_salmonella - densidad_ecoli

    comparacion = {
        "ecoli": {
            "longitud_genoma_pb": long_ecoli,
            "contenido_gc_porcentaje": round(gc_ecoli, 2),
            "total_genes_cds": total_cds_ecoli,
            "pb_codificante": pb_codificante_ecoli,
            "densidad_genica_porcentaje": round(densidad_ecoli, 2),
            "tamano_promedio_gen_pb": round(statistics.mean([g["longitud_pb"] for g in genes_ecoli]), 2)
        },
        "salmonella": {
            "longitud_genoma_pb": long_salmonella,
            "contenido_gc_porcentaje": round(gc_salmonella, 2),
            "total_genes_cds": total_cds_salmonella,
            "pb_codificante": pb_codificante_salmonella,
            "densidad_genica_porcentaje": round(densidad_salmonella, 2),
            "tamano_promedio_gen_pb": round(statistics.mean([g["longitud_pb"] for g in genes_salmonella]), 2)
        },
        "diferencias": {
            "longitud_pb": dif_longitud,
            "longitud_porcentaje": round((dif_longitud / long_ecoli) * 100, 2),
            "contenido_gc": round(dif_gc, 2),
            "total_genes": dif_genes,
            "densidad_genica": round(dif_densidad, 2)
        }
    }

    # Mostrar resultados
    print(f"\n  {'Metrica':<35} {'E. coli':>15} {'Salmonella':>15} {'Diferencia':>15}")
    print("  " + "-" * 80)

    print(f"  {'Longitud del genoma (pb)':<35} {long_ecoli:>15,} {long_salmonella:>15,} {dif_longitud:>+15,}")
    print(f"  {'Contenido GC (%)':<35} {gc_ecoli:>15.2f} {gc_salmonella:>15.2f} {dif_gc:>+15.2f}")
    print(f"  {'Total genes (CDS)':<35} {total_cds_ecoli:>15,} {total_cds_salmonella:>15,} {dif_genes:>+15,}")
    print(f"  {'Densidad genica (%)':<35} {densidad_ecoli:>15.2f} {densidad_salmonella:>15.2f} {dif_densidad:>+15.2f}")

    print(f"\n  [INTERPRETACION]")
    print(f"  - Salmonella tiene un genoma {abs(dif_longitud):,} pb ({abs(comparacion['diferencias']['longitud_porcentaje']):.1f}%) "
          f"{'mayor' if dif_longitud > 0 else 'menor'} que E. coli")
    print(f"  - Salmonella tiene {abs(dif_genes):,} genes {'mas' if dif_genes > 0 else 'menos'} que E. coli")
    print(f"  - El contenido GC es similar (diferencia de {abs(dif_gc):.2f}%)")

    return comparacion


def comparar_virulencia(genes_ecoli, genes_salmonella):
    """
    Compara genes de virulencia entre ambos organismos.

    Returns:
        dict: Comparacion de genes de virulencia
    """
    print("\n" + "=" * 70)
    print("COMPARACION DE GENES DE VIRULENCIA")
    print("=" * 70)
    print("  (Genes relacionados con patogenicidad, invasion, toxinas, etc.)")

    virulencia_ecoli = contar_genes_virulencia(genes_ecoli)
    virulencia_salmonella = contar_genes_virulencia(genes_salmonella)

    comparacion = {
        "ecoli": virulencia_ecoli,
        "salmonella": virulencia_salmonella,
        "diferencia_total": virulencia_salmonella["total"] - virulencia_ecoli["total"]
    }

    print(f"\n  CONTEO TOTAL DE GENES DE VIRULENCIA:")
    print(f"    E. coli K-12:     {virulencia_ecoli['total']:,} genes")
    print(f"    Salmonella LT2:   {virulencia_salmonella['total']:,} genes")
    print(f"    Diferencia:       {comparacion['diferencia_total']:+,} genes")

    print(f"\n  DESGLOSE POR CATEGORIA EN SALMONELLA:")
    for categoria, cantidad in sorted(virulencia_salmonella["categorias"].items(), key=lambda x: -x[1]):
        cant_ecoli = virulencia_ecoli["categorias"].get(categoria, 0)
        print(f"    {categoria:<35} {cantidad:>5} (E.coli: {cant_ecoli})")

    print(f"\n  [INTERPRETACION]")
    print(f"  - Salmonella tiene {comparacion['diferencia_total']:+} genes de virulencia mas que E. coli K-12")
    print(f"  - E. coli K-12 es una cepa de laboratorio NO patogena")
    print(f"  - Los genes de virulencia en E. coli K-12 son mayormente residuales/inactivos")
    print(f"  - Salmonella tiene genes activos para invadir celulas y causar enfermedad")

    if virulencia_salmonella["ejemplos"]:
        print(f"\n  EJEMPLOS DE GENES DE VIRULENCIA EN SALMONELLA (top 10):")
        for i, gen in enumerate(virulencia_salmonella["ejemplos"][:10], 1):
            print(f"    {i:2}. {gen['locus']}: {gen['producto'][:60]}")

    return comparacion


def comparar_distancias_intergenicas(genes_ecoli, genes_salmonella):
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

    dist_ecoli = calcular_distancias_intergenicas(genes_ecoli)
    dist_salmonella = calcular_distancias_intergenicas(genes_salmonella)

    comparacion = {
        "ecoli": dist_ecoli,
        "salmonella": dist_salmonella
    }

    print(f"\n  {'Metrica':<40} {'E. coli':>15} {'Salmonella':>15}")
    print("  " + "-" * 70)
    print(f"  {'Distancia intergenica promedio (pb)':<40} {dist_ecoli['promedio_pb']:>15.2f} {dist_salmonella['promedio_pb']:>15.2f}")
    print(f"  {'Distancia intergenica mediana (pb)':<40} {dist_ecoli['mediana_pb']:>15.2f} {dist_salmonella['mediana_pb']:>15.2f}")
    print(f"  {'Distancia maxima (pb)':<40} {dist_ecoli['maxima_pb']:>15,} {dist_salmonella['maxima_pb']:>15,}")
    print(f"  {'Regiones intergenicas > 5 kb':<40} {dist_ecoli['regiones_grandes_5kb']:>15} {dist_salmonella['regiones_grandes_5kb']:>15}")

    print(f"\n  [INTERPRETACION]")
    dif_regiones = dist_salmonella['regiones_grandes_5kb'] - dist_ecoli['regiones_grandes_5kb']
    if dif_regiones > 0:
        print(f"  - Salmonella tiene {dif_regiones} regiones intergenicas grandes mas que E. coli")
        print(f"  - Estas regiones pueden corresponder a ISLAS DE PATOGENICIDAD (SPIs)")
        print(f"  - Las SPIs contienen genes de virulencia adquiridos horizontalmente")
    else:
        print(f"  - Ambos genomas tienen distribucion similar de distancias intergenicas")

    if dist_salmonella['detalle_regiones_grandes']:
        print(f"\n  REGIONES INTERGENICAS GRANDES EN SALMONELLA (posibles SPIs):")
        for i, region in enumerate(dist_salmonella['detalle_regiones_grandes'][:5], 1):
            print(f"    {i}. Entre {region['gen_anterior']} y {region['gen_siguiente']}")
            print(f"       Posicion: {region['posicion_inicio']:,} - {region['posicion_fin']:,} ({region['distancia_pb']:,} pb)")

    return comparacion


def generar_resumen_comparativo(metricas, virulencia, distancias):
    """
    Genera un resumen narrativo de la comparacion.

    Returns:
        dict: Resumen con conclusiones principales
    """
    print("\n" + "=" * 70)
    print("RESUMEN: POR QUE SALMONELLA ES PATOGENA Y E. COLI K-12 NO")
    print("=" * 70)

    resumen = {
        "titulo": "Diferencias geneticas entre E. coli K-12 y Salmonella LT2",
        "divergencia_evolutiva": "~100-160 millones de anos",
        "similitud_genica": "~70% de genes compartidos",
        "diferencias_clave": []
    }

    # Diferencia 1: Tamano del genoma
    dif_tamano = metricas["diferencias"]["longitud_pb"]
    resumen["diferencias_clave"].append({
        "aspecto": "Tamano del genoma",
        "observacion": f"Salmonella tiene {abs(dif_tamano):,} pb adicionales",
        "significado": "ADN extra contiene genes de virulencia y adaptacion"
    })

    # Diferencia 2: Genes de virulencia
    dif_virulencia = virulencia["diferencia_total"]
    resumen["diferencias_clave"].append({
        "aspecto": "Genes de virulencia",
        "observacion": f"Salmonella tiene {dif_virulencia:+} genes de virulencia",
        "significado": "Genes para invadir celulas, secretar toxinas, evadir inmunidad"
    })

    # Diferencia 3: Islas de patogenicidad
    dif_islas = distancias["salmonella"]["regiones_grandes_5kb"] - distancias["ecoli"]["regiones_grandes_5kb"]
    resumen["diferencias_clave"].append({
        "aspecto": "Islas de patogenicidad (SPIs)",
        "observacion": f"Salmonella tiene {distancias['salmonella']['regiones_grandes_5kb']} regiones >5kb",
        "significado": "SPIs contienen sistemas de secrecion tipo III y efectores"
    })

    print(f"\n  CONTEXTO EVOLUTIVO:")
    print(f"    E. coli y Salmonella comparten un ancestro comun")
    print(f"    Divergieron hace {resumen['divergencia_evolutiva']}")
    print(f"    Comparten aproximadamente {resumen['similitud_genica']}")

    print(f"\n  DIFERENCIAS CLAVE QUE HACEN A SALMONELLA PATOGENA:")
    for i, dif in enumerate(resumen["diferencias_clave"], 1):
        print(f"\n    {i}. {dif['aspecto'].upper()}")
        print(f"       Observacion: {dif['observacion']}")
        print(f"       Significado: {dif['significado']}")

    print(f"\n  MECANISMOS DE PATOGENICIDAD DE SALMONELLA:")
    print(f"    1. INVASION: Sistema de secrecion tipo III (T3SS) inyecta proteinas")
    print(f"       en celulas del intestino para inducir su propia fagocitosis")
    print(f"    2. SUPERVIVENCIA: Sobrevive dentro de macrofagos evadiendo muerte")
    print(f"    3. DISEMINACION: Puede pasar del intestino al torrente sanguineo")
    print(f"    4. INFLAMACION: Induce respuesta inflamatoria que causa diarrea")

    print(f"\n  POR QUE E. COLI K-12 NO ES PATOGENA:")
    print(f"    1. Es una cepa de laboratorio derivada de un aislado fecal de 1922")
    print(f"    2. Ha perdido muchos genes de virulencia por domesticacion")
    print(f"    3. No tiene sistemas de secrecion tipo III funcionales")
    print(f"    4. No puede invadir celulas epiteliales")
    print(f"    5. Es la cepa modelo mas segura para investigacion")

    resumen["conclusiones"] = [
        "Pequenas diferencias geneticas (~5-10% del genoma) separan comensal de patogeno",
        "La adquisicion horizontal de islas de patogenicidad es clave en evolucion bacteriana",
        "E. coli K-12 y Salmonella son excelentes modelos para estudiar patogenicidad",
        "La comparacion genomica permite identificar genes de virulencia"
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


def exportar_comparacion_csv(metricas, nombre_archivo):
    """Exporta las metricas comparativas a CSV."""
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")

    with open(ruta, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Metrica', 'E_coli_K12', 'Salmonella_LT2', 'Diferencia', 'Unidad'])

        # Extraer metricas
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


def exportar_genes_virulencia_csv(virulencia_salmonella, nombre_archivo):
    """Exporta lista de genes de virulencia de Salmonella a CSV."""
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")

    with open(ruta, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Locus_tag', 'Nombre_gen', 'Producto'])

        for gen in virulencia_salmonella["ejemplos"]:
            writer.writerow([gen['locus'], gen['gen'], gen['producto']])

    print(f"       [OK] Guardado: {nombre_archivo}.csv")


# FUNCION PRINCIPAL
# =============================================================================

def main():
    """Funcion principal que ejecuta la comparacion de genomas."""

    print("\n" + "=" * 70)
    print("COMPARACION DE GENOMAS: E. coli vs Salmonella")
    print("=" * 70)
    print(f"\n  Fecha de analisis: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Objetivo: Comparar genomas para entender diferencias entre")
    print(f"            comensal (E. coli K-12) y patogeno (Salmonella LT2)")
    print("")

    # Crear directorios
    crear_directorios()

    # Verificar archivos
    archivos_ok, faltantes = verificar_archivos()
    if not archivos_ok:
        print("[ERROR] Faltan archivos GenBank:")
        for f in faltantes:
            print(f"        - {f}")
        print("\n[INFO] Ejecuta descargar_genoma.py y selecciona opcion 3 (ambos genomas)")
        return

    # Cargar genomas
    print("\n[PASO 1/4] Cargando archivos GenBank...")
    registro_ecoli = cargar_genbank("ecoli")
    registro_salmonella = cargar_genbank("salmonella")

    # Extraer genes
    print("\n[PASO 2/4] Extrayendo genes de ambos genomas...")
    genes_ecoli = extraer_genes_con_info(registro_ecoli, "ecoli")
    genes_salmonella = extraer_genes_con_info(registro_salmonella, "salmonella")
    print(f"       E. coli:     {len(genes_ecoli):,} genes CDS extraidos")
    print(f"       Salmonella:  {len(genes_salmonella):,} genes CDS extraidos")

    # Realizar comparaciones
    print("\n[PASO 3/4] Realizando comparaciones...")

    metricas = comparar_metricas_generales(
        registro_ecoli, registro_salmonella,
        genes_ecoli, genes_salmonella
    )

    virulencia = comparar_virulencia(genes_ecoli, genes_salmonella)

    distancias = comparar_distancias_intergenicas(genes_ecoli, genes_salmonella)

    resumen = generar_resumen_comparativo(metricas, virulencia, distancias)

    # Compilar resultados
    resultados_completos = {
        "fecha_analisis": datetime.now().isoformat(),
        "organismos_comparados": {
            "organismo_1": ORGANISMOS["ecoli"],
            "organismo_2": ORGANISMOS["salmonella"]
        },
        "metricas_generales": metricas,
        "genes_virulencia": virulencia,
        "distancias_intergenicas": distancias,
        "resumen_interpretativo": resumen
    }

    # Exportar resultados
    print("\n[PASO 4/4] Exportando resultados...")
    print("  (Los archivos se guardan en resultados/tablas/)\n")

    print("  [1/3] Exportando comparacion completa (JSON)...")
    exportar_comparacion_json(resultados_completos, "comparacion_ecoli_vs_salmonella")

    print("  [2/3] Exportando metricas comparativas (CSV)...")
    exportar_comparacion_csv(metricas, "metricas_comparacion")

    print("  [3/3] Exportando genes de virulencia Salmonella (CSV)...")
    exportar_genes_virulencia_csv(virulencia["salmonella"], "genes_virulencia_salmonella")

    # Resumen final
    print("\n" + "=" * 70)
    print("ARCHIVOS GENERADOS")
    print("=" * 70)
    print(f"    - comparacion_ecoli_vs_salmonella.json   (comparacion completa)")
    print(f"    - metricas_comparacion.csv              (metricas lado a lado)")
    print(f"    - genes_virulencia_salmonella.csv       (genes de patogenicidad)")

    print("\n" + "=" * 70)
    print("[OK] Comparacion de genomas completada exitosamente")
    print("=" * 70)
    print("\n  Siguiente paso: Ejecuta visualizaciones.py para generar graficos")
    print("                  de comparacion entre ambos genomas.\n")


if __name__ == "__main__":
    main()

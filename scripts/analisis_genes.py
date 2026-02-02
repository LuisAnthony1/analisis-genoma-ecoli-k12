#!/usr/bin/env python3
"""
Analisis de Genes del Genoma de E. coli K-12 MG1655

Este script extrae y analiza informacion de los genes anotados
en el archivo GenBank del genoma de E. coli K-12 MG1655.

Funcionalidades:
- Extraer todos los genes y CDS (secuencias codificantes) del GenBank
- Contar el numero total de genes (~4,300)
- Calcular la densidad genica del genoma
- Analizar la distribucion de tamanos de genes
- Calcular el contenido GC de las regiones codificantes
- Comparar resultados con valores de literatura cientifica
- Exportar resultados en formato CSV y JSON

Fecha: 2026
Proyecto: Analisis del genoma de E. coli K-12 MG1655
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

# Rutas de archivos
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_BASE, "datos", "crudo")
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")
ARCHIVO_GENBANK = os.path.join(RUTA_DATOS_CRUDO, "ecoli_k12.gb")

# Valores de referencia de literatura cientifica para E. coli K-12
VALORES_LITERATURA = {
    "longitud_genoma": 4641652,
    "genes_totales": 4300,          # Aproximadamente
    "genes_codificantes": 4290,     # CDS aproximados
    "contenido_gc_genoma": 50.8,    # Porcentaje
    "contenido_gc_cds": 51.5,       # GC en regiones codificantes
    "densidad_genica": 87.0,        # Porcentaje del genoma que codifica
    "tamano_promedio_gen": 950,     # pb aproximado
}


# FUNCIONES DE UTILIDAD
# =============================================================================

def crear_directorios():
    """Crea los directorios de resultados si no existen."""
    if not os.path.exists(RUTA_RESULTADOS):
        os.makedirs(RUTA_RESULTADOS)
        print(f"[INFO] Directorio creado: {RUTA_RESULTADOS}")


def cargar_genbank(ruta_archivo):
    """
    Carga el archivo GenBank del genoma.

    Args:
        ruta_archivo: Ruta al archivo GenBank

    Returns:
        SeqRecord: Objeto con el genoma y anotaciones
    """
    print(f"[INFO] Cargando archivo GenBank: {ruta_archivo}")

    if not os.path.exists(ruta_archivo):
        raise FileNotFoundError(f"No se encontro el archivo: {ruta_archivo}")

    registro = SeqIO.read(ruta_archivo, "genbank")

    print(f"[OK] GenBank cargado: {len(registro.seq):,} pb, {len(registro.features):,} features")
    return registro


# FUNCIONES DE EXTRACCION DE GENES
# =============================================================================

def extraer_genes(registro):
    """
    Extrae todos los genes del archivo GenBank.

    Args:
        registro: Objeto SeqRecord con el genoma

    Returns:
        list: Lista de diccionarios con informacion de cada gen
    """
    print("[INFO] Extrayendo genes del GenBank...")

    genes = []
    secuencia_genoma = registro.seq

    for feature in registro.features:
        if feature.type == "CDS":
            # Extraer informacion basica
            inicio = int(feature.location.start)
            fin = int(feature.location.end)
            hebra = feature.location.strand  # 1 = forward, -1 = reverse
            longitud = fin - inicio

            # Extraer secuencia del gen
            secuencia_gen = feature.location.extract(secuencia_genoma)

            # Calcular contenido GC del gen
            gc_gen = gc_fraction(secuencia_gen) * 100

            # Extraer anotaciones
            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
            nombre_gen = feature.qualifiers.get("gene", [""])[0]
            producto = feature.qualifiers.get("product", [""])[0]
            proteina_id = feature.qualifiers.get("protein_id", [""])[0]

            # Calcular numero de aminoacidos
            num_aminoacidos = longitud // 3

            gen_info = {
                "locus_tag": locus_tag,
                "nombre_gen": nombre_gen,
                "inicio": inicio,
                "fin": fin,
                "longitud_pb": longitud,
                "num_aminoacidos": num_aminoacidos,
                "hebra": "+" if hebra == 1 else "-",
                "contenido_gc": round(gc_gen, 2),
                "producto": producto,
                "proteina_id": proteina_id
            }

            genes.append(gen_info)

    print(f"[OK] Genes extraidos: {len(genes):,} CDS")
    return genes


def extraer_otros_features(registro):
    """
    Extrae y cuenta otros tipos de features del GenBank.

    Args:
        registro: Objeto SeqRecord con el genoma

    Returns:
        dict: Conteo de cada tipo de feature
    """
    conteo_features = Counter()

    for feature in registro.features:
        conteo_features[feature.type] += 1

    return dict(conteo_features)


# FUNCIONES DE ANALISIS
# =============================================================================

def analizar_estadisticas_genes(genes, longitud_genoma):
    """
    Calcula estadisticas generales de los genes.

    Args:
        genes: Lista de genes extraidos
        longitud_genoma: Longitud total del genoma

    Returns:
        dict: Estadisticas de los genes
    """
    print("\n" + "=" * 60)
    print("ESTADISTICAS GENERALES DE GENES")
    print("=" * 60)

    total_genes = len(genes)
    longitudes = [g["longitud_pb"] for g in genes]
    contenidos_gc = [g["contenido_gc"] for g in genes]

    # Calcular cobertura del genoma
    total_pb_codificante = sum(longitudes)
    densidad_genica = (total_pb_codificante / longitud_genoma) * 100

    # Estadisticas de tamano
    tamano_min = min(longitudes)
    tamano_max = max(longitudes)
    tamano_promedio = statistics.mean(longitudes)
    tamano_mediana = statistics.median(longitudes)
    desviacion_std = statistics.stdev(longitudes)

    # Estadisticas de GC
    gc_promedio = statistics.mean(contenidos_gc)
    gc_min = min(contenidos_gc)
    gc_max = max(contenidos_gc)

    # Genes por hebra
    genes_forward = sum(1 for g in genes if g["hebra"] == "+")
    genes_reverse = sum(1 for g in genes if g["hebra"] == "-")

    resultados = {
        "total_genes": total_genes,
        "longitud_genoma_pb": longitud_genoma,
        "total_pb_codificante": total_pb_codificante,
        "densidad_genica_porcentaje": round(densidad_genica, 2),
        "tamano_gen": {
            "minimo_pb": tamano_min,
            "maximo_pb": tamano_max,
            "promedio_pb": round(tamano_promedio, 2),
            "mediana_pb": round(tamano_mediana, 2),
            "desviacion_std": round(desviacion_std, 2)
        },
        "contenido_gc_cds": {
            "promedio": round(gc_promedio, 2),
            "minimo": round(gc_min, 2),
            "maximo": round(gc_max, 2)
        },
        "distribucion_hebras": {
            "forward": genes_forward,
            "reverse": genes_reverse,
            "porcentaje_forward": round((genes_forward / total_genes) * 100, 2),
            "porcentaje_reverse": round((genes_reverse / total_genes) * 100, 2)
        }
    }

    # Mostrar resultados
    print(f"\n  Total de genes (CDS):       {total_genes:,}")
    print(f"  Longitud del genoma:        {longitud_genoma:,} pb")
    print(f"  Total pb codificantes:      {total_pb_codificante:,} pb")
    print(f"  Densidad genica:            {densidad_genica:.2f}%")

    print(f"\n  TAMANO DE GENES:")
    print(f"    Minimo:     {tamano_min:,} pb")
    print(f"    Maximo:     {tamano_max:,} pb")
    print(f"    Promedio:   {tamano_promedio:,.2f} pb")
    print(f"    Mediana:    {tamano_mediana:,.2f} pb")
    print(f"    Desv. Std:  {desviacion_std:,.2f} pb")

    print(f"\n  CONTENIDO GC EN CDS:")
    print(f"    Promedio:   {gc_promedio:.2f}%")
    print(f"    Rango:      {gc_min:.2f}% - {gc_max:.2f}%")

    print(f"\n  DISTRIBUCION POR HEBRA:")
    print(f"    Forward (+): {genes_forward:,} ({(genes_forward/total_genes)*100:.1f}%)")
    print(f"    Reverse (-): {genes_reverse:,} ({(genes_reverse/total_genes)*100:.1f}%)")

    return resultados


def analizar_distribucion_tamanos(genes):
    """
    Analiza la distribucion de tamanos de genes por rangos.

    Args:
        genes: Lista de genes extraidos

    Returns:
        dict: Distribucion de tamanos por rangos
    """
    print("\n" + "=" * 60)
    print("DISTRIBUCION DE TAMANOS DE GENES")
    print("=" * 60)

    # Definir rangos de tamano
    rangos = [
        (0, 300, "0-300 pb (muy cortos)"),
        (301, 600, "301-600 pb (cortos)"),
        (601, 900, "601-900 pb (medianos)"),
        (901, 1500, "901-1500 pb (largos)"),
        (1501, 3000, "1501-3000 pb (muy largos)"),
        (3001, float('inf'), ">3000 pb (extra largos)")
    ]

    distribucion = {}
    total = len(genes)

    print(f"\n  {'Rango':<30} {'Cantidad':>10} {'Porcentaje':>12}")
    print("  " + "-" * 55)

    for min_val, max_val, etiqueta in rangos:
        cantidad = sum(1 for g in genes if min_val <= g["longitud_pb"] <= max_val)
        porcentaje = (cantidad / total) * 100
        distribucion[etiqueta] = {
            "cantidad": cantidad,
            "porcentaje": round(porcentaje, 2)
        }
        print(f"  {etiqueta:<30} {cantidad:>10,} {porcentaje:>11.2f}%")

    print("  " + "-" * 55)
    print(f"  {'TOTAL':<30} {total:>10,} {'100.00':>11}%")

    return distribucion


def comparar_con_literatura(estadisticas):
    """
    Compara los resultados con valores de literatura cientifica.

    Args:
        estadisticas: Diccionario con estadisticas calculadas

    Returns:
        dict: Comparacion con literatura
    """
    print("\n" + "=" * 60)
    print("COMPARACION CON LITERATURA CIENTIFICA")
    print("=" * 60)

    comparaciones = {}

    metricas = [
        ("Total genes", estadisticas["total_genes"], VALORES_LITERATURA["genes_totales"], "genes"),
        ("Densidad genica", estadisticas["densidad_genica_porcentaje"], VALORES_LITERATURA["densidad_genica"], "%"),
        ("GC en CDS", estadisticas["contenido_gc_cds"]["promedio"], VALORES_LITERATURA["contenido_gc_cds"], "%"),
        ("Tamano promedio", estadisticas["tamano_gen"]["promedio_pb"], VALORES_LITERATURA["tamano_promedio_gen"], "pb"),
    ]

    print(f"\n  {'Metrica':<20} {'Observado':>12} {'Literatura':>12} {'Diferencia':>12}")
    print("  " + "-" * 60)

    for nombre, observado, literatura, unidad in metricas:
        diferencia = observado - literatura
        diferencia_pct = (diferencia / literatura) * 100 if literatura != 0 else 0
        signo = "+" if diferencia > 0 else ""

        comparaciones[nombre] = {
            "observado": observado,
            "literatura": literatura,
            "diferencia": round(diferencia, 2),
            "diferencia_porcentaje": round(diferencia_pct, 2)
        }

        print(f"  {nombre:<20} {observado:>10.2f} {unidad} {literatura:>10.2f} {unidad} {signo}{diferencia:>10.2f}")

    print("\n  [INTERPRETACION]")

    # Evaluar diferencias
    dif_genes = abs(estadisticas["total_genes"] - VALORES_LITERATURA["genes_totales"])
    if dif_genes < 100:
        print("  [OK] El numero de genes coincide con la literatura (~4,300)")
    else:
        print(f"  [INFO] Diferencia de {dif_genes} genes con literatura")

    dif_densidad = abs(estadisticas["densidad_genica_porcentaje"] - VALORES_LITERATURA["densidad_genica"])
    if dif_densidad < 5:
        print("  [OK] La densidad genica es consistente (~87% del genoma codifica)")
    else:
        print(f"  [INFO] Diferencia de {dif_densidad:.1f}% en densidad genica")

    return comparaciones


# FUNCIONES DE EXPORTACION
# =============================================================================

def exportar_json(datos, nombre_archivo):
    """Exporta los resultados a formato JSON."""
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.json")

    with open(ruta, 'w', encoding='utf-8') as archivo:
        json.dump(datos, archivo, indent=2, ensure_ascii=False)

    print(f"[OK] Exportado: {ruta}")


def exportar_genes_csv(genes, nombre_archivo):
    """
    Exporta la lista de genes a formato CSV.

    Args:
        genes: Lista de genes
        nombre_archivo: Nombre del archivo (sin extension)
    """
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")

    with open(ruta, 'w', newline='', encoding='utf-8') as archivo:
        if genes:
            escritor = csv.DictWriter(archivo, fieldnames=genes[0].keys())
            escritor.writeheader()
            escritor.writerows(genes)

    print(f"[OK] Exportado: {ruta}")


def exportar_estadisticas_csv(estadisticas, comparaciones, nombre_archivo):
    """
    Exporta las estadisticas a formato CSV.

    Args:
        estadisticas: Diccionario con estadisticas
        comparaciones: Diccionario con comparaciones
        nombre_archivo: Nombre del archivo (sin extension)
    """
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")

    with open(ruta, 'w', newline='', encoding='utf-8') as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow(['Metrica', 'Valor', 'Unidad', 'Valor_Literatura', 'Diferencia_Porcentaje'])

        # Metricas principales
        filas = [
            ('Total_genes', estadisticas['total_genes'], 'genes',
             VALORES_LITERATURA['genes_totales'],
             comparaciones.get('Total genes', {}).get('diferencia_porcentaje', 0)),

            ('Densidad_genica', estadisticas['densidad_genica_porcentaje'], '%',
             VALORES_LITERATURA['densidad_genica'],
             comparaciones.get('Densidad genica', {}).get('diferencia_porcentaje', 0)),

            ('GC_promedio_CDS', estadisticas['contenido_gc_cds']['promedio'], '%',
             VALORES_LITERATURA['contenido_gc_cds'],
             comparaciones.get('GC en CDS', {}).get('diferencia_porcentaje', 0)),

            ('Tamano_promedio_gen', estadisticas['tamano_gen']['promedio_pb'], 'pb',
             VALORES_LITERATURA['tamano_promedio_gen'],
             comparaciones.get('Tamano promedio', {}).get('diferencia_porcentaje', 0)),

            ('Tamano_mediana_gen', estadisticas['tamano_gen']['mediana_pb'], 'pb', '-', '-'),
            ('Tamano_minimo_gen', estadisticas['tamano_gen']['minimo_pb'], 'pb', '-', '-'),
            ('Tamano_maximo_gen', estadisticas['tamano_gen']['maximo_pb'], 'pb', '-', '-'),

            ('Genes_forward', estadisticas['distribucion_hebras']['forward'], 'genes', '-', '-'),
            ('Genes_reverse', estadisticas['distribucion_hebras']['reverse'], 'genes', '-', '-'),

            ('Total_pb_codificante', estadisticas['total_pb_codificante'], 'pb', '-', '-'),
        ]

        for fila in filas:
            escritor.writerow(fila)

    print(f"[OK] Exportado: {ruta}")


# FUNCION PRINCIPAL
# =============================================================================

def main():
    """Funcion principal que ejecuta todo el analisis de genes."""

    print("\n" + "=" * 60)
    print("ANALISIS DE GENES - E. coli K-12 MG1655")
    print("=" * 60)
    print(f"Fecha de analisis: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Crear directorios
    crear_directorios()

    # Cargar GenBank
    try:
        registro = cargar_genbank(ARCHIVO_GENBANK)
    except FileNotFoundError as error:
        print(f"[ERROR] {error}")
        print("[INFO] Ejecuta primero descargar_genoma.py para obtener el genoma")
        return

    longitud_genoma = len(registro.seq)

    # Extraer genes
    genes = extraer_genes(registro)

    # Contar otros features
    otros_features = extraer_otros_features(registro)
    print(f"\n[INFO] Tipos de features en el GenBank:")
    for tipo, cantidad in sorted(otros_features.items(), key=lambda x: -x[1])[:10]:
        print(f"       {tipo}: {cantidad:,}")

    # Analizar estadisticas
    estadisticas = analizar_estadisticas_genes(genes, longitud_genoma)

    # Analizar distribucion de tamanos
    distribucion = analizar_distribucion_tamanos(genes)

    # Comparar con literatura
    comparaciones = comparar_con_literatura(estadisticas)

    # Compilar todos los resultados
    todos_resultados = {
        "fecha_analisis": datetime.now().isoformat(),
        "archivo_fuente": ARCHIVO_GENBANK,
        "estadisticas_generales": estadisticas,
        "distribucion_tamanos": distribucion,
        "comparacion_literatura": comparaciones,
        "tipos_features": otros_features,
        "valores_literatura": VALORES_LITERATURA
    }

    # Exportar resultados
    print("\n" + "=" * 60)
    print("EXPORTANDO RESULTADOS")
    print("=" * 60)

    exportar_json(todos_resultados, "analisis_genes_completo")
    exportar_genes_csv(genes, "lista_genes")
    exportar_estadisticas_csv(estadisticas, comparaciones, "estadisticas_genes")

    # Resumen final
    print("\n" + "=" * 60)
    print("RESUMEN FINAL")
    print("=" * 60)
    print(f"  Genoma analizado:       E. coli K-12 MG1655")
    print(f"  Total genes (CDS):      {estadisticas['total_genes']:,}")
    print(f"  Genes esperados:        ~{VALORES_LITERATURA['genes_totales']:,}")
    print(f"  Densidad genica:        {estadisticas['densidad_genica_porcentaje']:.2f}%")
    print(f"  GC en CDS:              {estadisticas['contenido_gc_cds']['promedio']:.2f}%")
    print(f"  Tamano promedio gen:    {estadisticas['tamano_gen']['promedio_pb']:,.0f} pb")
    print("=" * 60)

    print("\n[OK] Analisis completado exitosamente\n")


if __name__ == "__main__":
    main()

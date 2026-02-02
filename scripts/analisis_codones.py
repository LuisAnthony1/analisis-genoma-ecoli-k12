#!/usr/bin/env python3
"""
Analisis de Codones del Genoma de E. coli K-12 MG1655

Este script analiza los codones de inicio (ATG) y parada (TAA, TAG, TGA)
en el genoma completo de E. coli K-12 MG1655.

Funcionalidades:
- Contar todos los codones ATG y calcular densidad por kilobase
- Contar codones de parada (TAA, TAG, TGA) y sus proporciones
- Comparar resultados con valores de literatura cientifica
- Exportar resultados en formato CSV y JSON

Fecha: 2026
Proyecto: Analisis del genoma de E. coli K-12 MG1655
"""

from Bio import SeqIO
import os
import json
import csv
from datetime import datetime


# CONFIGURACION
# =============================================================================

# Rutas de archivos
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_BASE, "datos", "crudo")
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")
ARCHIVO_FASTA = os.path.join(RUTA_DATOS_CRUDO, "ecoli_k12.fasta")

# Codones a analizar
CODON_INICIO = "ATG"
CODONES_PARADA = ["TAA", "TAG", "TGA"]

# Valores de referencia de literatura cientifica para E. coli K-12
# Fuente: Blattner et al. (1997), GenBank, estudios de genomica comparativa
VALORES_LITERATURA = {
    "longitud_genoma": 4641652,
    "genes_anotados": 4300,  # Aproximadamente
    "contenido_gc": 50.8,   # Porcentaje
    # Proporciones de codones de parada en E. coli (literatura)
    # TAA es el mas comun (~63%), TAG (~8%), TGA (~29%)
    "proporcion_taa": 63.0,
    "proporcion_tag": 8.0,
    "proporcion_tga": 29.0,
}


# FUNCIONES DE UTILIDAD
# =============================================================================

def crear_directorios():
    """Crea los directorios de resultados si no existen."""
    if not os.path.exists(RUTA_RESULTADOS):
        os.makedirs(RUTA_RESULTADOS)
        print(f"[INFO] Directorio creado: {RUTA_RESULTADOS}")


def cargar_secuencia(ruta_archivo):
    """
    Carga la secuencia del genoma desde un archivo FASTA.

    Args:
        ruta_archivo: Ruta al archivo FASTA

    Returns:
        str: Secuencia de ADN en mayusculas
    """
    print(f"[INFO] Cargando secuencia desde: {ruta_archivo}")

    if not os.path.exists(ruta_archivo):
        raise FileNotFoundError(f"No se encontro el archivo: {ruta_archivo}")

    registro = SeqIO.read(ruta_archivo, "fasta")
    secuencia = str(registro.seq).upper()

    print(f"[OK] Secuencia cargada: {len(secuencia):,} pares de bases")
    return secuencia


# FUNCIONES DE ANALISIS
# =============================================================================

def contar_codon(secuencia, codon):
    """
    Cuenta las ocurrencias de un codon en la secuencia.

    Args:
        secuencia: Secuencia de ADN
        codon: Codon a buscar (3 nucleotidos)

    Returns:
        int: Numero de ocurrencias del codon
    """
    contador = 0
    codon = codon.upper()

    # Buscar el codon en toda la secuencia (no solo en marcos de lectura)
    for i in range(len(secuencia) - 2):
        if secuencia[i:i+3] == codon:
            contador += 1

    return contador


def calcular_densidad_por_kilobase(conteo, longitud_genoma):
    """
    Calcula la densidad de un codon por kilobase de genoma.

    Args:
        conteo: Numero de ocurrencias del codon
        longitud_genoma: Longitud total del genoma en pb

    Returns:
        float: Densidad por kilobase
    """
    kilobases = longitud_genoma / 1000
    return conteo / kilobases


def analizar_codones_inicio(secuencia):
    """
    Analiza los codones de inicio (ATG) en el genoma.

    Args:
        secuencia: Secuencia de ADN

    Returns:
        dict: Resultados del analisis de codones de inicio
    """
    print("\n" + "=" * 60)
    print("ANALISIS DE CODONES DE INICIO (ATG)")
    print("=" * 60)

    longitud = len(secuencia)
    conteo_atg = contar_codon(secuencia, CODON_INICIO)
    densidad = calcular_densidad_por_kilobase(conteo_atg, longitud)

    # Comparacion con genes anotados
    genes_anotados = VALORES_LITERATURA["genes_anotados"]
    ratio_atg_genes = conteo_atg / genes_anotados

    resultados = {
        "codon": CODON_INICIO,
        "conteo_total": conteo_atg,
        "longitud_genoma_pb": longitud,
        "densidad_por_kb": round(densidad, 4),
        "genes_anotados_literatura": genes_anotados,
        "ratio_atg_por_gen": round(ratio_atg_genes, 2),
    }

    # Mostrar resultados
    print(f"\n  Conteo total de ATG:        {conteo_atg:,}")
    print(f"  Longitud del genoma:        {longitud:,} pb")
    print(f"  Densidad por kilobase:      {densidad:.4f} ATG/kb")
    print(f"  Genes anotados (literatura): {genes_anotados:,}")
    print(f"  Ratio ATG/gen:              {ratio_atg_genes:.2f}")

    print("\n  [INTERPRETACION]")
    print(f"  Por cada gen anotado, hay aproximadamente {ratio_atg_genes:.0f} codones ATG")
    print("  en el genoma. Esto indica que NO todos los ATG corresponden")
    print("  a inicios reales de genes funcionales.")

    return resultados


def analizar_codones_parada(secuencia):
    """
    Analiza los codones de parada (TAA, TAG, TGA) en el genoma.

    Args:
        secuencia: Secuencia de ADN

    Returns:
        dict: Resultados del analisis de codones de parada
    """
    print("\n" + "=" * 60)
    print("ANALISIS DE CODONES DE PARADA (TAA, TAG, TGA)")
    print("=" * 60)

    longitud = len(secuencia)
    conteos = {}
    densidades = {}

    # Contar cada codon de parada
    for codon in CODONES_PARADA:
        conteos[codon] = contar_codon(secuencia, codon)
        densidades[codon] = calcular_densidad_por_kilobase(conteos[codon], longitud)

    # Calcular total y proporciones
    total_parada = sum(conteos.values())
    proporciones = {}
    for codon in CODONES_PARADA:
        proporciones[codon] = (conteos[codon] / total_parada) * 100

    # Construir resultados
    resultados = {
        "conteos": conteos,
        "total_codones_parada": total_parada,
        "densidades_por_kb": {k: round(v, 4) for k, v in densidades.items()},
        "proporciones_porcentaje": {k: round(v, 2) for k, v in proporciones.items()},
        "proporciones_literatura": {
            "TAA": VALORES_LITERATURA["proporcion_taa"],
            "TAG": VALORES_LITERATURA["proporcion_tag"],
            "TGA": VALORES_LITERATURA["proporcion_tga"],
        }
    }

    # Mostrar resultados
    print("\n  CONTEOS:")
    print(f"  {'Codon':<8} {'Conteo':>12} {'Densidad/kb':>15} {'Proporcion':>12}")
    print("  " + "-" * 50)
    for codon in CODONES_PARADA:
        print(f"  {codon:<8} {conteos[codon]:>12,} {densidades[codon]:>15.4f} {proporciones[codon]:>11.2f}%")
    print("  " + "-" * 50)
    print(f"  {'TOTAL':<8} {total_parada:>12,}")

    # Comparacion con literatura
    print("\n  COMPARACION CON LITERATURA CIENTIFICA:")
    print(f"  {'Codon':<8} {'Observado':>12} {'Literatura':>12} {'Diferencia':>12}")
    print("  " + "-" * 50)
    for codon in CODONES_PARADA:
        lit_key = f"proporcion_{codon.lower()}"
        valor_lit = VALORES_LITERATURA[lit_key]
        diferencia = proporciones[codon] - valor_lit
        signo = "+" if diferencia > 0 else ""
        print(f"  {codon:<8} {proporciones[codon]:>11.2f}% {valor_lit:>11.2f}% {signo}{diferencia:>11.2f}%")

    print("\n  [INTERPRETACION]")
    print("  En E. coli, TAA es el codon de parada preferido (~63% en genes),")
    print("  seguido de TGA (~29%) y TAG (~8%). Las proporciones observadas")
    print("  en el genoma completo difieren porque incluyen secuencias no codificantes.")

    return resultados


def calcular_contenido_gc(secuencia):
    """
    Calcula el contenido de GC del genoma.

    Args:
        secuencia: Secuencia de ADN

    Returns:
        dict: Resultados del contenido GC
    """
    print("\n" + "=" * 60)
    print("CONTENIDO DE GC")
    print("=" * 60)

    longitud = len(secuencia)
    conteo_g = secuencia.count('G')
    conteo_c = secuencia.count('C')
    conteo_a = secuencia.count('A')
    conteo_t = secuencia.count('T')

    contenido_gc = ((conteo_g + conteo_c) / longitud) * 100
    contenido_at = ((conteo_a + conteo_t) / longitud) * 100

    gc_literatura = VALORES_LITERATURA["contenido_gc"]
    diferencia = contenido_gc - gc_literatura

    resultados = {
        "conteo_nucleotidos": {
            "A": conteo_a,
            "T": conteo_t,
            "G": conteo_g,
            "C": conteo_c
        },
        "contenido_gc_porcentaje": round(contenido_gc, 2),
        "contenido_at_porcentaje": round(contenido_at, 2),
        "gc_literatura": gc_literatura,
        "diferencia_con_literatura": round(diferencia, 2)
    }

    print(f"\n  Nucleotidos:")
    print(f"    A: {conteo_a:,} ({(conteo_a/longitud)*100:.2f}%)")
    print(f"    T: {conteo_t:,} ({(conteo_t/longitud)*100:.2f}%)")
    print(f"    G: {conteo_g:,} ({(conteo_g/longitud)*100:.2f}%)")
    print(f"    C: {conteo_c:,} ({(conteo_c/longitud)*100:.2f}%)")
    print(f"\n  Contenido GC: {contenido_gc:.2f}%")
    print(f"  Contenido AT: {contenido_at:.2f}%")
    print(f"\n  Valor literatura: {gc_literatura}%")
    print(f"  Diferencia: {'+' if diferencia > 0 else ''}{diferencia:.2f}%")

    if abs(diferencia) < 1:
        print("\n  [OK] El contenido GC coincide con la literatura")
    else:
        print(f"\n  [ADVERTENCIA] Diferencia de {abs(diferencia):.2f}% con literatura")

    return resultados


# =============================================================================
# FUNCIONES DE EXPORTACION
# =============================================================================

def exportar_json(datos, nombre_archivo):
    """
    Exporta los resultados a formato JSON.

    Args:
        datos: Diccionario con los resultados
        nombre_archivo: Nombre del archivo (sin extension)
    """
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.json")

    with open(ruta, 'w', encoding='utf-8') as archivo:
        json.dump(datos, archivo, indent=2, ensure_ascii=False)

    print(f"[OK] Exportado: {ruta}")


def exportar_csv_codones(resultados_inicio, resultados_parada, nombre_archivo):
    """
    Exporta los resultados de codones a formato CSV.

    Args:
        resultados_inicio: Resultados del analisis de ATG
        resultados_parada: Resultados del analisis de codones de parada
        nombre_archivo: Nombre del archivo (sin extension)
    """
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")

    with open(ruta, 'w', newline='', encoding='utf-8') as archivo:
        escritor = csv.writer(archivo)

        # Encabezado
        escritor.writerow(['Tipo', 'Codon', 'Conteo', 'Densidad_por_kb', 'Proporcion_porcentaje'])

        # Codon de inicio
        escritor.writerow([
            'Inicio',
            resultados_inicio['codon'],
            resultados_inicio['conteo_total'],
            resultados_inicio['densidad_por_kb'],
            '100.00'
        ])

        # Codones de parada
        for codon in CODONES_PARADA:
            escritor.writerow([
                'Parada',
                codon,
                resultados_parada['conteos'][codon],
                resultados_parada['densidades_por_kb'][codon],
                resultados_parada['proporciones_porcentaje'][codon]
            ])

    print(f"[OK] Exportado: {ruta}")


def exportar_resumen_csv(todos_resultados, nombre_archivo):
    """
    Exporta un resumen general a CSV.

    Args:
        todos_resultados: Diccionario con todos los resultados
        nombre_archivo: Nombre del archivo (sin extension)
    """
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")

    with open(ruta, 'w', newline='', encoding='utf-8') as archivo:
        escritor = csv.writer(archivo)

        # Encabezado
        escritor.writerow(['Metrica', 'Valor', 'Unidad', 'Valor_Literatura', 'Diferencia'])

        # Metricas generales
        escritor.writerow([
            'Longitud_genoma',
            todos_resultados['codones_inicio']['longitud_genoma_pb'],
            'pb',
            VALORES_LITERATURA['longitud_genoma'],
            todos_resultados['codones_inicio']['longitud_genoma_pb'] - VALORES_LITERATURA['longitud_genoma']
        ])

        escritor.writerow([
            'Contenido_GC',
            todos_resultados['contenido_gc']['contenido_gc_porcentaje'],
            '%',
            VALORES_LITERATURA['contenido_gc'],
            todos_resultados['contenido_gc']['diferencia_con_literatura']
        ])

        escritor.writerow([
            'Total_ATG',
            todos_resultados['codones_inicio']['conteo_total'],
            'codones',
            '-',
            '-'
        ])

        escritor.writerow([
            'Densidad_ATG',
            todos_resultados['codones_inicio']['densidad_por_kb'],
            'ATG/kb',
            '-',
            '-'
        ])

        # Codones de parada
        for codon in CODONES_PARADA:
            lit_key = f"proporcion_{codon.lower()}"
            valor_obs = todos_resultados['codones_parada']['proporciones_porcentaje'][codon]
            valor_lit = VALORES_LITERATURA[lit_key]
            escritor.writerow([
                f'Proporcion_{codon}',
                valor_obs,
                '%',
                valor_lit,
                round(valor_obs - valor_lit, 2)
            ])

    print(f"[OK] Exportado: {ruta}")


# FUNCION PRINCIPAL
# =============================================================================

def main():
    """Funcion principal que ejecuta todo el analisis de codones."""

    print("\n" + "=" * 60)
    print("ANALISIS DE CODONES - E. coli K-12 MG1655")
    print("=" * 60)
    print(f"Fecha de analisis: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Crear directorios
    crear_directorios()

    # Cargar secuencia
    try:
        secuencia = cargar_secuencia(ARCHIVO_FASTA)
    except FileNotFoundError as error:
        print(f"[ERROR] {error}")
        print("[INFO] Ejecuta primero descargar_genoma.py para obtener el genoma")
        return

    # Realizar analisis
    resultados_inicio = analizar_codones_inicio(secuencia)
    resultados_parada = analizar_codones_parada(secuencia)
    resultados_gc = calcular_contenido_gc(secuencia)

    # Compilar todos los resultados
    todos_resultados = {
        "fecha_analisis": datetime.now().isoformat(),
        "archivo_fuente": ARCHIVO_FASTA,
        "codones_inicio": resultados_inicio,
        "codones_parada": resultados_parada,
        "contenido_gc": resultados_gc,
        "valores_literatura": VALORES_LITERATURA
    }

    # Exportar resultados
    print("\n" + "=" * 60)
    print("EXPORTANDO RESULTADOS")
    print("=" * 60)

    exportar_json(todos_resultados, "analisis_codones_completo")
    exportar_csv_codones(resultados_inicio, resultados_parada, "codones_conteo")
    exportar_resumen_csv(todos_resultados, "resumen_analisis")

    # Resumen final
    print("\n" + "=" * 60)
    print("RESUMEN FINAL")
    print("=" * 60)
    print(f"  Genoma analizado:    E. coli K-12 MG1655")
    print(f"  Longitud:            {len(secuencia):,} pb")
    print(f"  Contenido GC:        {resultados_gc['contenido_gc_porcentaje']}%")
    print(f"  Total ATG:           {resultados_inicio['conteo_total']:,}")
    print(f"  Densidad ATG:        {resultados_inicio['densidad_por_kb']:.4f} /kb")
    print(f"  Total codones parada: {resultados_parada['total_codones_parada']:,}")
    print("=" * 60)

    print("\n[OK] Analisis completado exitosamente\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Analisis de Genes de Genomas Bacterianos

Este script extrae y analiza informacion de los genes anotados
en archivos GenBank de genomas bacterianos.

Organismos soportados:
1. Escherichia coli K-12 MG1655 (cepa de laboratorio)
2. Salmonella enterica serovar Typhimurium LT2 (patogeno)

Funcionalidades:
- Extraer todos los genes y CDS (secuencias codificantes) del GenBank
- Contar el numero total de genes
- Calcular la densidad genica del genoma
- Analizar la distribucion de tamanos de genes
- Calcular el contenido GC de las regiones codificantes
- Comparar resultados con valores de literatura cientifica
- Exportar resultados en formato CSV y JSON

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

# Rutas de archivos
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_BASE, "datos", "crudo")
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")

# Configuracion de organismos disponibles
ORGANISMOS = {
    1: {
        "nombre": "Escherichia coli K-12 MG1655",
        "nombre_corto": "ecoli_k12",
        "archivo_genbank": "ecoli_k12.gb",
        "descripcion": "Cepa de laboratorio de referencia, no patogena",
        "valores_literatura": {
            "longitud_genoma": 4641652,       # pb - NCBI NC_000913.3
            "genes_totales": 4319,            # CDS segun EcoCyc 2024
            "genes_codificantes": 4319,       # CDS (Coding DNA Sequences) - EcoCyc
            "genes_totales_anotados": 4651,   # Incluye ARN y pseudogenes
            "contenido_gc_genoma": 50.79,     # Porcentaje - NCBI
            "contenido_gc_cds": 51.06,        # GC en regiones codificantes
            "densidad_genica": 87.8,          # Porcentaje del genoma que codifica
            "tamano_promedio_gen": 940,       # pb aproximado - EcoCyc
        },
        "referencias": {
            "ncbi_refseq": "NCBI Reference Sequence NC_000913.3",
            "database": "EcoCyc Database (https://ecocyc.org/)",
            "paper": "Blattner FR et al. (1997) Science 277:1453-1462",
        }
    },
    2: {
        "nombre": "Salmonella enterica serovar Typhimurium LT2",
        "nombre_corto": "salmonella_lt2",
        "archivo_genbank": "salmonella_lt2.gb",
        "descripcion": "Patogeno - causa gastroenteritis, cepa de referencia",
        "valores_literatura": {
            "longitud_genoma": 4857450,       # pb - NCBI NC_003197.2
            "genes_totales": 4525,            # CDS aproximados
            "genes_codificantes": 4525,       # CDS (Coding DNA Sequences)
            "genes_totales_anotados": 4747,   # Incluye ARN y pseudogenes
            "contenido_gc_genoma": 52.22,     # Porcentaje - NCBI
            "contenido_gc_cds": 52.5,         # GC en regiones codificantes
            "densidad_genica": 86.8,          # Porcentaje del genoma que codifica
            "tamano_promedio_gen": 945,       # pb aproximado
        },
        "referencias": {
            "ncbi_refseq": "NCBI Reference Sequence NC_003197.2",
            "database": "SalmonellaDB / KEGG",
            "paper": "McClelland M et al. (2001) Nature 413:852-856",
        }
    }
}

# Variables globales que se configuran segun el organismo seleccionado
ORGANISMO_ACTUAL = None
ARCHIVO_GENBANK = None
VALORES_LITERATURA = None
REFERENCIAS = None


def seleccionar_organismo():
    """
    Selecciona el organismo a analizar desde sys.argv[1] (basename del genoma).

    Returns:
        int o str: Numero del organismo si es conocido, o basename directamente
    """
    if len(sys.argv) < 2:
        print("[ERROR] Uso: python analisis_genes.py <genome_basename>")
        print("  Ejemplo: python analisis_genes.py ecoli_k12")
        sys.exit(1)

    basename = sys.argv[1]

    # Buscar en organismos conocidos por nombre_corto
    for num, org in ORGANISMOS.items():
        if org["nombre_corto"] == basename:
            return num

    # Genoma desconocido - retornar basename como string
    return basename


def configurar_organismo(seleccion):
    """
    Configura las variables globales segun el organismo seleccionado.

    Args:
        seleccion: Numero del organismo (1 o 2) o basename como string
    """
    global ORGANISMO_ACTUAL, ARCHIVO_GENBANK, VALORES_LITERATURA, REFERENCIAS

    if isinstance(seleccion, int) and seleccion in ORGANISMOS:
        # Organismo conocido
        config = ORGANISMOS[seleccion]
        ORGANISMO_ACTUAL = config
        ARCHIVO_GENBANK = os.path.join(RUTA_DATOS_CRUDO, config["archivo_genbank"])
        VALORES_LITERATURA = config["valores_literatura"]
        REFERENCIAS = config["referencias"]
        print(f"\n[OK] Organismo configurado: {config['nombre']}")
        print(f"     Archivo: {config['archivo_genbank']}")
    else:
        # Genoma desconocido - configurar con valores genericos
        basename = str(seleccion)
        ARCHIVO_GENBANK = os.path.join(RUTA_DATOS_CRUDO, f"{basename}.gb")
        VALORES_LITERATURA = {}
        REFERENCIAS = {"fuente": "Archivo GenBank local"}
        ORGANISMO_ACTUAL = {
            "nombre": basename.replace("_", " ").title(),
            "nombre_corto": basename,
            "archivo_genbank": f"{basename}.gb",
            "descripcion": "Genoma cargado desde archivo local",
            "valores_literatura": VALORES_LITERATURA,
            "referencias": REFERENCIAS
        }
        print(f"\n[OK] Genoma configurado: {basename}")
        print(f"     Archivo: {basename}.gb")


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

    print(f"[OK] GenBank cargado exitosamente:")
    print(f"     - Longitud del genoma: {len(registro.seq):,} pares de bases")
    print(f"     - Elementos anotados:  {len(registro.features):,} features (genes, ARN, etc.)")
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
                "proteina_id": proteina_id,
                "secuencia": str(secuencia_gen)  # Guardamos la secuencia completa
            }

            genes.append(gen_info)

    print(f"[OK] Extraccion completada: {len(genes):,} secuencias codificantes (CDS) encontradas")
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


def analizar_genes_extremos(genes):
    """
    Analiza y muestra informacion detallada del gen mas largo y mas corto.
    Incluye secuencia (primeros y ultimos 10 nucleotidos), proteina y posiciones.

    Args:
        genes: Lista de genes extraidos
    """
    print("\n" + "=" * 70)
    print("ANALISIS DE GENES EXTREMOS (MAS LARGO Y MAS CORTO)")
    print("=" * 70)
    print("  (Mostramos los genes con tamanos extremos para entender la variabilidad)")

    # Ordenar genes por longitud
    genes_ordenados = sorted(genes, key=lambda g: g["longitud_pb"], reverse=True)

    # Gen mas largo
    gen_largo = genes_ordenados[0]
    # Gen mas corto
    gen_corto = genes_ordenados[-1]

    # ========== GEN MAS LARGO ==========
    print("\n" + "-" * 70)
    print("  GEN MAS LARGO DEL GENOMA")
    print("-" * 70)

    secuencia_largo = gen_largo["secuencia"]
    primeros_10_largo = secuencia_largo[:10]
    ultimos_10_largo = secuencia_largo[-10:]

    print(f"\n  IDENTIFICACION:")
    print(f"    Locus tag:          {gen_largo['locus_tag']}")
    print(f"    Nombre del gen:     {gen_largo['nombre_gen'] if gen_largo['nombre_gen'] else '(sin nombre asignado)'}")
    print(f"    ID de proteina:     {gen_largo['proteina_id']}")

    print(f"\n  UBICACION EN EL GENOMA:")
    print(f"    Posicion de inicio: {gen_largo['inicio']:,} (par de bases donde comienza)")
    print(f"    Posicion de fin:    {gen_largo['fin']:,} (par de bases donde termina)")
    print(f"    Hebra:              {gen_largo['hebra']} ({'forward 5->3' if gen_largo['hebra'] == '+' else 'reverse 3->5'})")

    print(f"\n  TAMANO:")
    print(f"    Longitud total:     {gen_largo['longitud_pb']:,} pares de bases")
    print(f"    Aminoacidos:        {gen_largo['num_aminoacidos']:,} aminoacidos en la proteina")
    print(f"    Contenido GC:       {gen_largo['contenido_gc']:.2f}%")

    print(f"\n  SECUENCIA DE NUCLEOTIDOS:")
    print(f"    Primeros 10 nucleotidos: 5'-{primeros_10_largo}-... (inicio del gen)")
    print(f"    Ultimos 10 nucleotidos:  ...-{ultimos_10_largo}-3' (final del gen)")
    print(f"    Secuencia completa:      {len(secuencia_largo):,} nucleotidos totales")

    print(f"\n  PROTEINA QUE CODIFICA:")
    print(f"    Producto: {gen_largo['producto']}")
    if gen_largo['producto']:
        print(f"    (Esta proteina tiene {gen_largo['num_aminoacidos']:,} aminoacidos de longitud)")

    # ========== GEN MAS CORTO ==========
    print("\n" + "-" * 70)
    print("  GEN MAS CORTO DEL GENOMA")
    print("-" * 70)

    secuencia_corto = gen_corto["secuencia"]
    primeros_10_corto = secuencia_corto[:10] if len(secuencia_corto) >= 10 else secuencia_corto
    ultimos_10_corto = secuencia_corto[-10:] if len(secuencia_corto) >= 10 else secuencia_corto

    print(f"\n  IDENTIFICACION:")
    print(f"    Locus tag:          {gen_corto['locus_tag']}")
    print(f"    Nombre del gen:     {gen_corto['nombre_gen'] if gen_corto['nombre_gen'] else '(sin nombre asignado)'}")
    print(f"    ID de proteina:     {gen_corto['proteina_id']}")

    print(f"\n  UBICACION EN EL GENOMA:")
    print(f"    Posicion de inicio: {gen_corto['inicio']:,} (par de bases donde comienza)")
    print(f"    Posicion de fin:    {gen_corto['fin']:,} (par de bases donde termina)")
    print(f"    Hebra:              {gen_corto['hebra']} ({'forward 5->3' if gen_corto['hebra'] == '+' else 'reverse 3->5'})")

    print(f"\n  TAMANO:")
    print(f"    Longitud total:     {gen_corto['longitud_pb']:,} pares de bases")
    print(f"    Aminoacidos:        {gen_corto['num_aminoacidos']:,} aminoacidos en la proteina")
    print(f"    Contenido GC:       {gen_corto['contenido_gc']:.2f}%")

    print(f"\n  SECUENCIA DE NUCLEOTIDOS:")
    if len(secuencia_corto) <= 20:
        print(f"    Secuencia completa:      5'-{secuencia_corto}-3'")
        print(f"    (El gen es tan corto que mostramos la secuencia completa)")
    else:
        print(f"    Primeros 10 nucleotidos: 5'-{primeros_10_corto}-... (inicio del gen)")
        print(f"    Ultimos 10 nucleotidos:  ...-{ultimos_10_corto}-3' (final del gen)")
    print(f"    Total nucleotidos:       {len(secuencia_corto):,} nucleotidos")

    print(f"\n  PROTEINA QUE CODIFICA:")
    print(f"    Producto: {gen_corto['producto']}")
    if gen_corto['producto']:
        print(f"    (Esta proteina tiene solo {gen_corto['num_aminoacidos']:,} aminoacidos - muy pequena)")

    return {
        "gen_mas_largo": {
            "locus_tag": gen_largo["locus_tag"],
            "nombre": gen_largo["nombre_gen"],
            "longitud_pb": gen_largo["longitud_pb"],
            "inicio": gen_largo["inicio"],
            "fin": gen_largo["fin"],
            "producto": gen_largo["producto"],
            "primeros_10_nt": primeros_10_largo,
            "ultimos_10_nt": ultimos_10_largo
        },
        "gen_mas_corto": {
            "locus_tag": gen_corto["locus_tag"],
            "nombre": gen_corto["nombre_gen"],
            "longitud_pb": gen_corto["longitud_pb"],
            "inicio": gen_corto["inicio"],
            "fin": gen_corto["fin"],
            "producto": gen_corto["producto"],
            "primeros_10_nt": primeros_10_corto,
            "ultimos_10_nt": ultimos_10_corto
        }
    }


def analizar_genes_vs_cds(registro, genes_cds):
    """
    Analiza la diferencia entre todos los genes anotados y los CDS (codificantes).
    Explica cuantos genes fueron "descartados" y por que.

    Args:
        registro: Objeto SeqRecord con el genoma
        genes_cds: Lista de genes CDS extraidos

    Returns:
        dict: Analisis detallado de genes vs CDS
    """
    print("\n" + "=" * 70)
    print("ANALISIS: GENES TOTALES vs GENES CODIFICANTES (CDS)")
    print("=" * 70)
    print("  (Diferencia entre todos los genes y los que codifican proteinas)")

    # Contar todos los features tipo "gene"
    genes_totales = []
    genes_tipo_gene = []
    for feature in registro.features:
        if feature.type == "gene":
            locus = feature.qualifiers.get("locus_tag", [""])[0]
            nombre = feature.qualifiers.get("gene", [""])[0]
            genes_tipo_gene.append({
                "locus_tag": locus,
                "nombre": nombre,
                "inicio": int(feature.location.start),
                "fin": int(feature.location.end)
            })

    total_genes_anotados = len(genes_tipo_gene)
    total_cds = len(genes_cds)

    # Crear set de locus_tag de CDS para comparar
    locus_cds = set(g["locus_tag"] for g in genes_cds)

    # Encontrar genes que NO son CDS
    genes_no_codificantes = []
    for gene in genes_tipo_gene:
        if gene["locus_tag"] not in locus_cds:
            genes_no_codificantes.append(gene)

    # Contar otros tipos de RNA
    conteo_rna = Counter()
    genes_rna_info = []
    for feature in registro.features:
        if feature.type in ["tRNA", "rRNA", "ncRNA", "tmRNA"]:
            locus = feature.qualifiers.get("locus_tag", [""])[0]
            producto = feature.qualifiers.get("product", [""])[0]
            conteo_rna[feature.type] += 1
            genes_rna_info.append({
                "tipo": feature.type,
                "locus_tag": locus,
                "producto": producto
            })

    # Pseudogenes
    pseudogenes = []
    for feature in registro.features:
        if feature.type == "gene":
            if "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers:
                locus = feature.qualifiers.get("locus_tag", [""])[0]
                pseudogenes.append(locus)

    print(f"\n  CONTEO GENERAL DE GENES:")
    print(f"  " + "-" * 65)
    print(f"    Total de genes anotados en GenBank:     {total_genes_anotados:,} genes")
    print(f"    Genes que codifican proteinas (CDS):   {total_cds:,} genes")
    print(f"    Diferencia (genes no codificantes):    {total_genes_anotados - total_cds:,} genes")

    print(f"\n  DESGLOSE DE GENES NO CODIFICANTES:")
    print(f"  " + "-" * 65)
    print(f"    (Estos genes NO producen proteinas, pero tienen otras funciones)")
    print(f"")

    for tipo_rna, cantidad in sorted(conteo_rna.items(), key=lambda x: -x[1]):
        descripcion = {
            "tRNA": "ARN de transferencia - transportan aminoacidos durante sintesis de proteinas",
            "rRNA": "ARN ribosomal - componentes estructurales de los ribosomas",
            "ncRNA": "ARN no codificante - regulan expresion genica",
            "tmRNA": "ARN de transferencia-mensajero - rescata ribosomas atascados"
        }.get(tipo_rna, "Otro tipo de ARN")
        print(f"    {tipo_rna}: {cantidad:,} genes")
        print(f"         Funcion: {descripcion}")
        print()

    if pseudogenes:
        print(f"    Pseudogenes: {len(pseudogenes):,} genes")
        print(f"         Funcion: Genes inactivos/no funcionales (mutaciones acumuladas)")
        print()

    print(f"\n  RESUMEN EXPLICATIVO:")
    print(f"  " + "-" * 65)
    print(f"    [1] GENES TOTALES ({total_genes_anotados:,}):")
    print(f"        Todos los segmentos de ADN identificados como 'genes' en el GenBank.")
    print(f"        Incluye tanto genes codificantes como no codificantes.")
    print(f"")
    print(f"    [2] GENES CODIFICANTES - CDS ({total_cds:,}):")
    print(f"        Solo los genes que producen proteinas (Coding DNA Sequences).")
    print(f"        Segun EcoCyc (2024): {VALORES_LITERATURA['genes_codificantes']:,} CDS en E. coli K-12 MG1655.")
    print(f"        Cada uno tiene: codon de inicio (ATG), codones, codon de parada.")
    print(f"")
    print(f"    [3] GENES NO CODIFICANTES ({total_genes_anotados - total_cds:,}):")
    print(f"        Genes que producen ARN funcional pero NO proteinas:")
    total_rna = sum(conteo_rna.values())
    print(f"        - ARN funcionales (tRNA, rRNA, etc.): {total_rna:,} genes")
    if pseudogenes:
        print(f"        - Pseudogenes (genes inactivos): {len(pseudogenes):,} genes")
    print(f"")
    print(f"    [IMPORTANTE] Los {total_cds:,} CDS son los genes 'utiles' para producir")
    print(f"    proteinas. Los otros {total_genes_anotados - total_cds:,} genes tienen funciones regulatorias")
    print(f"    o estructurales (ARN) pero no se 'desechan' - son igualmente importantes.")

    return {
        "genes_totales_anotados": total_genes_anotados,
        "genes_codificantes_cds": total_cds,
        "genes_no_codificantes": total_genes_anotados - total_cds,
        "desglose_rna": dict(conteo_rna),
        "pseudogenes": len(pseudogenes),
        "explicacion": "Los genes no codificantes producen ARN funcional (tRNA, rRNA) necesarios para la vida celular"
    }


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
    print("=" * 70)

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
    print(f"\n  Total de genes (CDS):       {total_genes:,} secuencias codificantes de proteinas")
    print(f"  Longitud del genoma:        {longitud_genoma:,} pares de bases (nucleotidos)")
    print(f"  Total pares de bases que codifican proteinas: {total_pb_codificante:,} pares de bases")
    print(f"  Densidad genica:            {densidad_genica:.2f}% del genoma codifica proteinas")

    print(f"\n  TAMANO DE GENES (en pares de bases - cada 3 pares = 1 aminoacido):")
    print(f"    Gen mas corto:      {tamano_min:,} pares de bases ({tamano_min//3} aminoacidos)")
    print(f"    Gen mas largo:      {tamano_max:,} pares de bases ({tamano_max//3} aminoacidos)")
    print(f"    Tamano promedio:    {tamano_promedio:,.2f} pares de bases ({tamano_promedio/3:.0f} aminoacidos)")
    print(f"    Tamano mediana:     {tamano_mediana:,.2f} pares de bases (valor central)")
    print(f"    Desviacion estandar: {desviacion_std:,.2f} pares de bases (variabilidad)")

    print(f"\n  CONTENIDO GC EN REGIONES CODIFICANTES (CDS):")
    print(f"    (GC = porcentaje de Guanina + Citosina en el ADN)")
    print(f"    Promedio:   {gc_promedio:.2f}% de G+C")
    print(f"    Rango:      {gc_min:.2f}% (minimo) - {gc_max:.2f}% (maximo)")

    print(f"\n  DISTRIBUCION POR HEBRA DEL ADN:")
    print(f"    (El ADN tiene dos hebras: forward 5'->3' y reverse 3'->5')")
    print(f"    Hebra Forward (+): {genes_forward:,} genes ({(genes_forward/total_genes)*100:.1f}%)")
    print(f"    Hebra Reverse (-): {genes_reverse:,} genes ({(genes_reverse/total_genes)*100:.1f}%)")

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
    print("=" * 70)
    print("  (Clasificacion de genes segun su longitud en pares de bases)")
    print("  (Recordar: 3 pares de bases = 1 codon = 1 aminoacido)")

    # Definir rangos de tamano con descripciones detalladas
    rangos = [
        (0, 300, "0-300 pares de bases", "muy cortos, <100 aminoacidos"),
        (301, 600, "301-600 pares de bases", "cortos, 100-200 aminoacidos"),
        (601, 900, "601-900 pares de bases", "medianos, 200-300 aminoacidos"),
        (901, 1500, "901-1500 pares de bases", "largos, 300-500 aminoacidos"),
        (1501, 3000, "1501-3000 pares de bases", "muy largos, 500-1000 aminoacidos"),
        (3001, float('inf'), ">3000 pares de bases", "extra largos, >1000 aminoacidos")
    ]

    distribucion = {}
    total = len(genes)

    print(f"\n  {'Rango de tamano':<28} {'Descripcion':<30} {'Genes':>8} {'%':>8}")
    print("  " + "-" * 78)

    for min_val, max_val, etiqueta, descripcion in rangos:
        cantidad = sum(1 for g in genes if min_val <= g["longitud_pb"] <= max_val)
        porcentaje = (cantidad / total) * 100
        distribucion[etiqueta] = {
            "cantidad": cantidad,
            "porcentaje": round(porcentaje, 2),
            "descripcion": descripcion
        }
        print(f"  {etiqueta:<28} {descripcion:<30} {cantidad:>8,} {porcentaje:>7.2f}%")

    print("  " + "-" * 78)
    print(f"  {'TOTAL DE GENES':<58} {total:>8,} {'100.00':>7}%")

    return distribucion


def comparar_con_literatura(estadisticas):
    """
    Compara los resultados con valores de literatura cientifica.

    Args:
        estadisticas: Diccionario con estadisticas calculadas

    Returns:
        dict: Comparacion con literatura
    """
    print("\n" + "=" * 70)
    print("COMPARACION CON LITERATURA CIENTIFICA")
    print("=" * 70)

    if not VALORES_LITERATURA:
        print("  No hay valores de literatura disponibles para este organismo.")
        return {}

    print("  (Validacion de resultados comparando con valores publicados)")
    print("")
    print("  FUENTES DE REFERENCIA:")
    if REFERENCIAS:
        for clave, valor in REFERENCIAS.items():
            print(f"    - {clave}: {valor}")

    comparaciones = {}

    metricas = [
        ("Total CDS", estadisticas["total_genes"], VALORES_LITERATURA.get("genes_codificantes"), "genes", "EcoCyc"),
        ("Densidad genica", estadisticas["densidad_genica_porcentaje"], VALORES_LITERATURA.get("densidad_genica"), "%", "EcoCyc"),
        ("GC en CDS", estadisticas["contenido_gc_cds"]["promedio"], VALORES_LITERATURA.get("contenido_gc_cds"), "%", "NCBI"),
        ("Tamano promedio", estadisticas["tamano_gen"]["promedio_pb"], VALORES_LITERATURA.get("tamano_promedio_gen"), "pb", "EcoCyc"),
    ]

    # Filtrar metricas sin valores de literatura
    metricas = [(n, o, l, u, f) for n, o, l, u, f in metricas if l is not None]

    print(f"\n  {'Metrica':<18} {'Analisis':>12} {'Literatura':>12} {'Dif.':>10} {'Fuente':<10}")
    print("  " + "-" * 70)

    for nombre, observado, literatura, unidad, fuente in metricas:
        diferencia = observado - literatura
        diferencia_pct = (diferencia / literatura) * 100 if literatura != 0 else 0
        signo = "+" if diferencia > 0 else ""

        comparaciones[nombre] = {
            "observado": observado,
            "literatura": literatura,
            "diferencia": round(diferencia, 2),
            "diferencia_porcentaje": round(diferencia_pct, 2),
            "fuente": fuente
        }

        print(f"  {nombre:<18} {observado:>10.2f} {unidad:<3} {literatura:>10.2f} {unidad:<3} {signo}{diferencia:>8.2f}  [{fuente}]")

    print("\n  [INTERPRETACION DE RESULTADOS]")

    # Evaluar diferencias
    dif_genes = abs(estadisticas["total_genes"] - VALORES_LITERATURA["genes_codificantes"])
    if dif_genes < 10:
        print(f"  [OK] EXCELENTE: El numero de CDS ({estadisticas['total_genes']:,}) coincide con EcoCyc ({VALORES_LITERATURA['genes_codificantes']:,})")
    elif dif_genes < 50:
        print(f"  [OK] El numero de genes es muy cercano a la literatura (diferencia: {dif_genes})")
    else:
        print(f"  [INFO] Diferencia de {dif_genes} genes con literatura - verificar version del GenBank")

    dif_densidad = abs(estadisticas["densidad_genica_porcentaje"] - VALORES_LITERATURA["densidad_genica"])
    if dif_densidad < 2:
        print(f"  [OK] La densidad genica ({estadisticas['densidad_genica_porcentaje']:.2f}%) es consistente con literatura ({VALORES_LITERATURA['densidad_genica']}%)")
    else:
        print(f"  [INFO] Diferencia de {dif_densidad:.1f}% en densidad genica")

    print("")
    print("  [NOTA] Pequenas variaciones son normales debido a:")
    print("         - Actualizaciones en la anotacion del genoma (NCBI actualiza periodicamente)")
    print("         - Diferentes criterios para clasificar genes/pseudogenes")
    print("         - Version del archivo GenBank descargado")

    return comparaciones


# FUNCIONES DE EXPORTACION
# =============================================================================

def exportar_json(datos, nombre_archivo):
    """Exporta los resultados a formato JSON."""
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.json")

    with open(ruta, 'w', encoding='utf-8') as archivo:
        json.dump(datos, archivo, indent=2, ensure_ascii=False)

    print(f"       [OK] Guardado: {nombre_archivo}.json")
    print(f"           Contiene: estadisticas, distribuciones y comparaciones")


def exportar_genes_csv(genes, nombre_archivo):
    """
    Exporta la lista de genes a formato CSV.
    Nota: Se excluye la secuencia completa para reducir el tamano del archivo.

    Args:
        genes: Lista de genes
        nombre_archivo: Nombre del archivo (sin extension)
    """
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")

    # Excluir la secuencia completa del CSV (muy grande)
    campos = [k for k in genes[0].keys() if k != "secuencia"]

    with open(ruta, 'w', newline='', encoding='utf-8') as archivo:
        if genes:
            escritor = csv.DictWriter(archivo, fieldnames=campos, extrasaction='ignore')
            escritor.writeheader()
            escritor.writerows(genes)

    print(f"       [OK] Guardado: {nombre_archivo}.csv")
    print(f"           Contiene: {len(genes):,} genes con locus_tag, nombre, posicion, tamano, GC, producto")


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

    print(f"       [OK] Guardado: {nombre_archivo}.csv")
    print(f"           Contiene: metricas principales comparadas con literatura cientifica")


# FUNCION PRINCIPAL
# =============================================================================

def main():
    """Funcion principal que ejecuta todo el analisis de genes."""

    # Paso 0: Seleccionar organismo desde sys.argv
    seleccion = seleccionar_organismo()
    configurar_organismo(seleccion)

    nombre_organismo = ORGANISMO_ACTUAL["nombre"]
    nombre_corto = ORGANISMO_ACTUAL["nombre_corto"]

    print("\n" + "=" * 70)
    print(f"ANALISIS DE GENES DEL GENOMA DE {nombre_organismo}")
    print("=" * 70)
    print(f"\n  Fecha de analisis: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Objetivo: Extraer y analizar todos los genes codificantes del genoma")
    print(f"  Metodo: Lectura del archivo GenBank y extraccion de CDS (coding sequences)")
    print("")

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
    print(f"\n[INFO] Tipos de elementos anotados en el archivo GenBank:")
    print(f"       (Features = elementos biologicos identificados en el genoma)")
    for tipo, cantidad in sorted(otros_features.items(), key=lambda x: -x[1])[:10]:
        descripcion_tipo = {
            "CDS": "secuencias codificantes de proteinas",
            "gene": "genes identificados",
            "source": "fuente del genoma",
            "rRNA": "ARN ribosomal",
            "tRNA": "ARN de transferencia",
            "misc_feature": "otras caracteristicas",
            "repeat_region": "regiones repetitivas",
            "mobile_element": "elementos moviles (transposones)"
        }.get(tipo, "")
        if descripcion_tipo:
            print(f"       {tipo}: {cantidad:,} ({descripcion_tipo})")
        else:
            print(f"       {tipo}: {cantidad:,}")

    # Verificar que se encontraron genes
    if len(genes) == 0:
        print("\n" + "=" * 70)
        print("[ERROR] NO SE ENCONTRARON GENES CDS EN EL ARCHIVO GENBANK")
        print("=" * 70)
        print("  El archivo GenBank parece no tener anotaciones de genes.")
        print("  Esto puede ocurrir si:")
        print("    1. El archivo se descargo sin anotaciones (solo secuencia)")
        print("    2. El formato de descarga no fue el correcto")
        print("")
        print("  SOLUCION: Vuelve a ejecutar descargar_genoma.py")
        print("            y selecciona el organismo nuevamente.")
        print("=" * 70 + "\n")
        return

    # Analizar genes totales vs CDS (codificantes vs no codificantes)
    analisis_genes_cds = analizar_genes_vs_cds(registro, genes)

    # Analizar genes extremos (mas largo y mas corto)
    genes_extremos = analizar_genes_extremos(genes)

    # Analizar estadisticas
    estadisticas = analizar_estadisticas_genes(genes, longitud_genoma)

    # Analizar distribucion de tamanos
    distribucion = analizar_distribucion_tamanos(genes)

    # Comparar con literatura
    comparaciones = comparar_con_literatura(estadisticas)

    # Compilar todos los resultados (sin incluir secuencias completas para reducir tamano)
    todos_resultados = {
        "fecha_analisis": datetime.now().isoformat(),
        "organismo": nombre_organismo,
        "organismo_corto": nombre_corto,
        "archivo_fuente": ARCHIVO_GENBANK,
        "estadisticas_generales": estadisticas,
        "distribucion_tamanos": distribucion,
        "comparacion_literatura": comparaciones,
        "tipos_features": otros_features,
        "valores_literatura": VALORES_LITERATURA,
        "referencias_bibliograficas": REFERENCIAS,
        "analisis_genes_vs_cds": analisis_genes_cds,
        "genes_extremos": genes_extremos,
        "nota": "Las secuencias completas de genes no se incluyen para reducir el tamano del archivo"
    }

    # Exportar resultados (con nombre del organismo en el archivo)
    print("\n" + "=" * 60)
    print("EXPORTANDO RESULTADOS A ARCHIVOS")
    print("=" * 70)
    print("  (Los archivos se guardan en la carpeta resultados/tablas/)\n")

    # Nombres de archivos con prefijo del organismo
    archivo_json = f"analisis_genes_{nombre_corto}"
    archivo_lista = f"lista_genes_{nombre_corto}"
    archivo_stats = f"estadisticas_genes_{nombre_corto}"

    print("  [1/3] Exportando analisis completo en formato JSON...")
    exportar_json(todos_resultados, archivo_json)

    print("  [2/3] Exportando lista de genes en formato CSV...")
    exportar_genes_csv(genes, archivo_lista)

    print("  [3/3] Exportando estadisticas en formato CSV...")
    exportar_estadisticas_csv(estadisticas, comparaciones, archivo_stats)

    # Resumen final
    print("\n" + "=" * 70)
    print("RESUMEN FINAL DEL ANALISIS")
    print("=" * 70)
    print(f"\n  ORGANISMO ANALIZADO:")
    print(f"    Nombre:       {nombre_organismo}")
    print(f"    Descripcion:  {ORGANISMO_ACTUAL['descripcion']}")

    print(f"\n  CONTEO DE GENES:")
    print(f"    Genes totales anotados:       {analisis_genes_cds['genes_totales_anotados']:,} genes en el GenBank")
    print(f"    Genes codificantes (CDS):     {analisis_genes_cds['genes_codificantes_cds']:,} genes que producen proteinas")
    print(f"    Genes no codificantes:        {analisis_genes_cds['genes_no_codificantes']:,} genes (ARN funcional)")
    print(f"    CDS segun literatura:         {VALORES_LITERATURA['genes_codificantes']:,} CDS")

    print(f"\n  GENES EXTREMOS:")
    print(f"    Gen mas largo:  {genes_extremos['gen_mas_largo']['locus_tag']} ({genes_extremos['gen_mas_largo']['longitud_pb']:,} pares de bases)")
    print(f"                    Proteina: {genes_extremos['gen_mas_largo']['producto'][:50]}...")
    print(f"    Gen mas corto:  {genes_extremos['gen_mas_corto']['locus_tag']} ({genes_extremos['gen_mas_corto']['longitud_pb']:,} pares de bases)")
    print(f"                    Proteina: {genes_extremos['gen_mas_corto']['producto']}")

    print(f"\n  METRICAS PRINCIPALES:")
    print(f"    Densidad genica:              {estadisticas['densidad_genica_porcentaje']:.2f}% del genoma codifica proteinas")
    print(f"    Contenido GC en genes:        {estadisticas['contenido_gc_cds']['promedio']:.2f}% de Guanina+Citosina")
    print(f"    Tamano promedio de genes:     {estadisticas['tamano_gen']['promedio_pb']:,.0f} pares de bases ({estadisticas['tamano_gen']['promedio_pb']/3:.0f} aminoacidos)")

    print(f"\n  ARCHIVOS GENERADOS:")
    print(f"    - {archivo_lista}.csv:         Lista completa de {estadisticas['total_genes']:,} genes")
    print(f"    - {archivo_stats}.csv:    Metricas resumidas")
    print(f"    - {archivo_json}.json: Todos los datos")

    print("\n" + "=" * 70)
    print(f"[OK] Analisis de {nombre_corto} completado exitosamente")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()

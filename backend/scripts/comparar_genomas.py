#!/usr/bin/env python3
"""
Comparacion Avanzada de Genomas Bacterianos

Script completo para comparar dos genomas bacterianos a partir de sus
archivos GenBank. Incluye analisis de ortologos, mutaciones, sintenia,
islas genomicas, perfil funcional, tRNA/rRNA, resistencia antibiotica,
ANI estimado y GC regional.

Uso: python comparar_genomas.py <basename_1> <basename_2>
"""

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
from Bio.Data.CodonTable import standard_dna_table
import os
import sys
import json
import csv
from datetime import datetime
from collections import Counter, defaultdict
import statistics
import random
import math

# CONFIGURACION
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIRECTORIO_PROYECTO = os.path.dirname(DIRECTORIO_BASE)
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_PROYECTO, "datos", "crudo")
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")

if len(sys.argv) < 3:
    print("[ERROR] Uso: python comparar_genomas.py <basename_1> <basename_2>")
    sys.exit(1)

BASENAME_1 = sys.argv[1]
BASENAME_2 = sys.argv[2]
ARCHIVO_GB_1 = os.path.join(RUTA_DATOS_CRUDO, f"{BASENAME_1}.gb")
ARCHIVO_GB_2 = os.path.join(RUTA_DATOS_CRUDO, f"{BASENAME_2}.gb")

PALABRAS_VIRULENCIA = [
    "invasion", "invasin", "virulence", "pathogen", "toxin", "secretion",
    "type III", "T3SS", "SPI", "effector", "adhesin", "fimbri", "pili",
    "flagell", "motility", "iron", "siderophore", "enterobactin",
    "resistance", "antibiotic", "drug"
]

# Tabla de codones para traduccion
CODON_TABLE = standard_dna_table.forward_table
STOP_CODONS = standard_dna_table.stop_codons

# Categorias funcionales por palabras clave en el producto
CATEGORIAS_FUNCIONALES = {
    "Metabolismo": ["dehydrogenase", "kinase", "synthase", "synthetase", "reductase",
                     "oxidase", "transferase", "hydrolase", "isomerase", "ligase",
                     "carboxylase", "phosphatase", "mutase", "aldolase", "epimerase",
                     "decarboxylase", "aminotransferase", "transaldolase", "enolase"],
    "Transporte": ["transporter", "permease", "porin", "channel", "pump", "export",
                    "import", "ABC", "symporter", "antiporter", "carrier", "efflux"],
    "Regulacion": ["regulator", "repressor", "activator", "transcription", "sigma",
                    "response regulator", "sensor", "two-component", "anti-sigma"],
    "Replicacion/Reparacion": ["DNA polymerase", "helicase", "primase", "ligase",
                                "topoisomerase", "gyrase", "recombinase", "repair",
                                "endonuclease", "exonuclease", "SSB", "recA"],
    "Traduccion/Ribosoma": ["ribosomal protein", "tRNA", "aminoacyl", "translation",
                             "elongation factor", "initiation factor", "release factor"],
    "Pared celular/Membrana": ["murein", "peptidoglycan", "lipopolysaccharide", "LPS",
                                "outer membrane", "inner membrane", "lipoprotein",
                                "penicillin-binding", "transpeptidase"],
    "Movilidad/Quimiotaxis": ["flagell", "motility", "chemotaxis", "che", "fli", "flg",
                               "motor", "hook", "basal body", "swimming"],
    "Defensa/Resistencia": ["resistance", "antibiotic", "drug", "beta-lactamase",
                             "efflux", "multidrug", "toxin-antitoxin"],
    "Virulencia/Patogenicidad": ["virulence", "pathogen", "invasion", "toxin",
                                  "secretion system", "effector", "adhesin"],
    "Estres/Chaperones": ["chaperone", "heat shock", "cold shock", "stress",
                           "GroEL", "GroES", "DnaK", "DnaJ", "ClpB", "IbpA"],
    "Proteinas hipoteticas": ["hypothetical", "uncharacterized", "DUF", "unknown function"],
}

# Genes de resistencia antibiotica (ampliado)
GENES_RESISTENCIA = {
    "Beta-lactamicos": ["beta-lactamase", "penicillin", "cephalosporin", "carbapenem",
                         "ampC", "bla", "OXA", "CTX-M", "TEM", "SHV", "KPC", "NDM"],
    "Aminoglucosidos": ["aminoglycoside", "streptomycin", "kanamycin", "gentamicin",
                         "aac", "aph", "ant", "rmt"],
    "Tetraciclinas": ["tetracycline", "tet(", "tetA", "tetB", "tetC", "tetR"],
    "Quinolonas": ["quinolone", "gyrase", "topoisomerase IV", "qnr", "fluoroquinolone"],
    "Macrolidos": ["macrolide", "erythromycin", "erm", "mef", "mph"],
    "Sulfonamidas/Trimetoprim": ["sulfonamide", "dihydrofolate", "trimethoprim", "sul1", "sul2", "dfr"],
    "Cloranfenicol": ["chloramphenicol", "cat", "florfenicol", "cfr"],
    "Polimixinas": ["polymyxin", "colistin", "mcr-", "arnT", "pmrA", "pmrB"],
    "Multidrogas": ["multidrug", "efflux pump", "mar", "sox", "acrAB", "tolC"],
    "Vancomicina": ["vancomycin", "vanA", "vanB", "vanC"],
}


# =============================================================================
# FUNCIONES DE CARGA
# =============================================================================

def crear_directorios():
    if not os.path.exists(RUTA_RESULTADOS):
        os.makedirs(RUTA_RESULTADOS)

def verificar_archivos():
    faltantes = []
    if not os.path.exists(ARCHIVO_GB_1):
        faltantes.append(f"{BASENAME_1}.gb")
    if not os.path.exists(ARCHIVO_GB_2):
        faltantes.append(f"{BASENAME_2}.gb")
    return len(faltantes) == 0, faltantes

def cargar_genbank(archivo_gb, nombre):
    print(f"[INFO] Cargando {nombre}...")
    registro = SeqIO.read(archivo_gb, "genbank")
    print(f"       Longitud: {len(registro.seq):,} pares de bases")
    return registro

def obtener_nombre_organismo(registro, basename):
    nombre = registro.annotations.get("organism", "")
    if nombre:
        return nombre
    return basename.replace("_", " ").title()


# =============================================================================
# FUNCIONES DE EXTRACCION
# =============================================================================

def extraer_genes_con_info(registro, nombre_organismo):
    """Extrae genes CDS con informacion detallada."""
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
            hebra = 1 if feature.location.strand == 1 else -1
            # Extraer secuencia de ADN
            seq_adn = str(secuencia_gen)
            es_virulencia = any(p.lower() in producto.lower() for p in PALABRAS_VIRULENCIA)
            genes.append({
                "locus_tag": locus_tag,
                "nombre_gen": nombre_gen,
                "inicio": inicio,
                "fin": fin,
                "longitud_pb": longitud,
                "contenido_gc": round(gc_gen, 2),
                "producto": producto,
                "hebra": hebra,
                "secuencia_adn": seq_adn,
                "es_virulencia": es_virulencia,
                "organismo": nombre_organismo
            })
    return genes


def extraer_features_no_cds(registro):
    """Extrae tRNA, rRNA y otros features no-CDS."""
    features = {"tRNA": [], "rRNA": [], "ncRNA": [], "misc_RNA": []}
    for feature in registro.features:
        ftype = feature.type
        if ftype in features:
            info = {
                "tipo": ftype,
                "inicio": int(feature.location.start),
                "fin": int(feature.location.end),
                "longitud": int(feature.location.end) - int(feature.location.start),
                "producto": feature.qualifiers.get("product", [""])[0],
                "nombre": feature.qualifiers.get("gene", [""])[0],
                "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
            }
            if ftype == "tRNA":
                info["anticodon"] = feature.qualifiers.get("anticodon", [""])[0]
                # Aminoacido del tRNA
                prod = info["producto"]
                if "tRNA-" in prod:
                    info["aminoacido"] = prod.replace("tRNA-", "").strip()
                else:
                    info["aminoacido"] = ""
            if ftype == "rRNA":
                info["subtipo"] = ""
                prod = info["producto"].lower()
                if "16s" in prod:
                    info["subtipo"] = "16S"
                elif "23s" in prod:
                    info["subtipo"] = "23S"
                elif "5s" in prod:
                    info["subtipo"] = "5S"
            features[ftype].append(info)
    return features


def calcular_distancias_intergenicas(genes):
    """Calcula distancias entre genes consecutivos."""
    genes_ordenados = sorted(genes, key=lambda g: g["inicio"])
    distancias = []
    regiones_grandes = []
    for i in range(1, len(genes_ordenados)):
        gen_anterior = genes_ordenados[i - 1]
        gen_actual = genes_ordenados[i]
        distancia = gen_actual["inicio"] - gen_anterior["fin"]
        if distancia > 0:
            distancias.append(distancia)
            if distancia > 5000:
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
        "detalle_regiones_grandes": regiones_grandes[:10]
    }


def contar_genes_virulencia(genes):
    """Cuenta genes relacionados con virulencia."""
    genes_virulencia = [g for g in genes if g["es_virulencia"]]
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


# =============================================================================
# COMPARACIONES BASICAS (existentes)
# =============================================================================

def comparar_metricas_generales(registro_1, registro_2, genes_1, genes_2, nombre_1, nombre_2):
    """Compara metricas generales entre ambos genomas."""
    print("\n" + "=" * 70)
    print("COMPARACION DE METRICAS GENERALES")
    print("=" * 70)

    long_1 = len(registro_1.seq)
    gc_1 = gc_fraction(registro_1.seq) * 100
    total_cds_1 = len(genes_1)
    pb_codificante_1 = sum(g["longitud_pb"] for g in genes_1)
    densidad_1 = (pb_codificante_1 / long_1) * 100 if long_1 > 0 else 0

    long_2 = len(registro_2.seq)
    gc_2 = gc_fraction(registro_2.seq) * 100
    total_cds_2 = len(genes_2)
    pb_codificante_2 = sum(g["longitud_pb"] for g in genes_2)
    densidad_2 = (pb_codificante_2 / long_2) * 100 if long_2 > 0 else 0

    dif_longitud = long_2 - long_1
    dif_gc = gc_2 - gc_1
    dif_genes = total_cds_2 - total_cds_1
    dif_densidad = densidad_2 - densidad_1

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

    n1_s = nombre_1[:20]
    n2_s = nombre_2[:20]
    print(f"\n  {'Metrica':<35} {n1_s:>20} {n2_s:>20} {'Diferencia':>15}")
    print("  " + "-" * 90)
    print(f"  {'Longitud (pb)':<35} {long_1:>20,} {long_2:>20,} {dif_longitud:>+15,}")
    print(f"  {'GC (%)':<35} {gc_1:>20.2f} {gc_2:>20.2f} {dif_gc:>+15.2f}")
    print(f"  {'Genes CDS':<35} {total_cds_1:>20,} {total_cds_2:>20,} {dif_genes:>+15,}")
    print(f"  {'Densidad (%)':<35} {densidad_1:>20.2f} {densidad_2:>20.2f} {dif_densidad:>+15.2f}")
    return comparacion


def comparar_virulencia(genes_1, genes_2, nombre_1, nombre_2):
    """Compara genes de virulencia."""
    print("\n" + "=" * 70)
    print("COMPARACION DE GENES DE VIRULENCIA")
    print("=" * 70)
    virulencia_1 = contar_genes_virulencia(genes_1)
    virulencia_2 = contar_genes_virulencia(genes_2)
    return {
        "ecoli": virulencia_1,
        "salmonella": virulencia_2,
        "diferencia_total": virulencia_2["total"] - virulencia_1["total"]
    }


def comparar_distancias_intergenicas(genes_1, genes_2, nombre_1, nombre_2):
    """Compara distancias intergenicas."""
    print("\n" + "=" * 70)
    print("COMPARACION DE DISTANCIAS INTERGENICAS")
    print("=" * 70)
    dist_1 = calcular_distancias_intergenicas(genes_1)
    dist_2 = calcular_distancias_intergenicas(genes_2)
    return {"ecoli": dist_1, "salmonella": dist_2}


def comparar_uso_codones(registro_1, registro_2, genes_1, genes_2, nombre_1, nombre_2):
    """Compara uso de codones."""
    print("\n" + "=" * 70)
    print("COMPARACION DE USO DE CODONES")
    print("=" * 70)

    def contar_codones(registro, genes):
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
    top10_1 = codones_1.most_common(10)
    top10_2 = codones_2.most_common(10)
    freq_1 = {c: round(n / total_1 * 100, 3) if total_1 > 0 else 0 for c, n in codones_1.items()}
    freq_2 = {c: round(n / total_2 * 100, 3) if total_2 > 0 else 0 for c, n in codones_2.items()}
    todos_codones = set(list(freq_1.keys()) + list(freq_2.keys()))
    diferencias = []
    for codon in todos_codones:
        f1 = freq_1.get(codon, 0)
        f2 = freq_2.get(codon, 0)
        diferencias.append({"codon": codon, "frecuencia_1": f1, "frecuencia_2": f2,
                            "diferencia_abs": round(abs(f2 - f1), 3)})
    diferencias.sort(key=lambda x: x["diferencia_abs"], reverse=True)
    return {
        "ecoli": {"total_codones": total_1,
                   "top_10": [{"codon": c, "conteo": n, "frecuencia": freq_1.get(c, 0)} for c, n in top10_1]},
        "salmonella": {"total_codones": total_2,
                        "top_10": [{"codon": c, "conteo": n, "frecuencia": freq_2.get(c, 0)} for c, n in top10_2]},
        "mayores_diferencias": diferencias[:15]
    }


def comparar_distribucion_tamanos(genes_1, genes_2, nombre_1, nombre_2):
    """Compara distribucion de tamanos de genes."""
    print("\n" + "=" * 70)
    print("COMPARACION DE DISTRIBUCION DE TAMANOS")
    print("=" * 70)

    def calc_stats(genes):
        tamanos = [g["longitud_pb"] for g in genes]
        if not tamanos:
            return {}
        return {
            "promedio": round(statistics.mean(tamanos), 1),
            "mediana": round(statistics.median(tamanos), 1),
            "minimo": min(tamanos), "maximo": max(tamanos),
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

    return {
        "ecoli": {"estadisticas": calc_stats(genes_1), "rangos": calc_rangos(genes_1)},
        "salmonella": {"estadisticas": calc_stats(genes_2), "rangos": calc_rangos(genes_2)}
    }


def comparar_distribucion_gc(genes_1, genes_2, nombre_1, nombre_2):
    """Compara distribucion de GC por gen."""
    print("\n" + "=" * 70)
    print("COMPARACION DE DISTRIBUCION GC")
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

    return {
        "ecoli": {"estadisticas": calc_gc_stats(genes_1), "rangos": calc_rangos_gc(genes_1)},
        "salmonella": {"estadisticas": calc_gc_stats(genes_2), "rangos": calc_rangos_gc(genes_2)}
    }


# =============================================================================
# NUEVOS ANALISIS AVANZADOS
# =============================================================================

def normalizar_producto(producto):
    """Normaliza el nombre de un producto para matching de ortologos."""
    p = producto.lower().strip()
    # Remover prefijos comunes que no afectan identidad
    for prefix in ["putative ", "probable ", "predicted ", "uncharacterized "]:
        p = p.replace(prefix, "")
    # Remover sufijos de organismo
    for suffix in [" protein", " family protein"]:
        if p.endswith(suffix) and len(p) > len(suffix) + 3:
            pass  # mantener - es parte del nombre
    return p.strip()


def detectar_ortologos(genes_1, genes_2, nombre_1, nombre_2):
    """
    Detecta genes ortologos (compartidos) y unicos entre dos genomas.
    Usa matching por nombre de producto normalizado.
    """
    print("\n" + "=" * 70)
    print("DETECCION DE ORTOLOGOS (GENES COMPARTIDOS vs UNICOS)")
    print("=" * 70)

    # Crear indices por producto normalizado
    idx_1 = defaultdict(list)
    idx_2 = defaultdict(list)

    for g in genes_1:
        if g["producto"] and "hypothetical" not in g["producto"].lower():
            key = normalizar_producto(g["producto"])
            if key:
                idx_1[key].append(g)

    for g in genes_2:
        if g["producto"] and "hypothetical" not in g["producto"].lower():
            key = normalizar_producto(g["producto"])
            if key:
                idx_2[key].append(g)

    # Tambien intentar match por nombre de gen
    genname_1 = {}
    genname_2 = {}
    for g in genes_1:
        if g["nombre_gen"]:
            genname_1[g["nombre_gen"].lower()] = g
    for g in genes_2:
        if g["nombre_gen"]:
            genname_2[g["nombre_gen"].lower()] = g

    # Encontrar compartidos
    productos_compartidos = set(idx_1.keys()) & set(idx_2.keys())
    genes_nombre_compartidos = set(genname_1.keys()) & set(genname_2.keys())

    # Construir pares de ortologos
    ortologos = []
    productos_ya_emparejados = set()

    for prod in productos_compartidos:
        g1_list = idx_1[prod]
        g2_list = idx_2[prod]
        # Emparejar el primer gen de cada lista
        for g1 in g1_list[:1]:
            for g2 in g2_list[:1]:
                ortologos.append({
                    "producto": g1["producto"],
                    "gen_1": {"locus": g1["locus_tag"], "nombre": g1["nombre_gen"],
                              "inicio": g1["inicio"], "fin": g1["fin"],
                              "longitud": g1["longitud_pb"], "gc": g1["contenido_gc"],
                              "hebra": g1["hebra"]},
                    "gen_2": {"locus": g2["locus_tag"], "nombre": g2["nombre_gen"],
                              "inicio": g2["inicio"], "fin": g2["fin"],
                              "longitud": g2["longitud_pb"], "gc": g2["contenido_gc"],
                              "hebra": g2["hebra"]},
                    "dif_longitud": abs(g1["longitud_pb"] - g2["longitud_pb"]),
                    "dif_gc": round(abs(g1["contenido_gc"] - g2["contenido_gc"]), 2),
                    "misma_hebra": g1["hebra"] == g2["hebra"],
                })
                productos_ya_emparejados.add(prod)

    # Agregar ortologos por nombre de gen que no fueron emparejados por producto
    for gname in genes_nombre_compartidos:
        g1 = genname_1[gname]
        g2 = genname_2[gname]
        prod_key = normalizar_producto(g1["producto"])
        if prod_key not in productos_ya_emparejados:
            ortologos.append({
                "producto": g1["producto"],
                "gen_1": {"locus": g1["locus_tag"], "nombre": g1["nombre_gen"],
                          "inicio": g1["inicio"], "fin": g1["fin"],
                          "longitud": g1["longitud_pb"], "gc": g1["contenido_gc"],
                          "hebra": g1["hebra"]},
                "gen_2": {"locus": g2["locus_tag"], "nombre": g2["nombre_gen"],
                          "inicio": g2["inicio"], "fin": g2["fin"],
                          "longitud": g2["longitud_pb"], "gc": g2["contenido_gc"],
                          "hebra": g2["hebra"]},
                "dif_longitud": abs(g1["longitud_pb"] - g2["longitud_pb"]),
                "dif_gc": round(abs(g1["contenido_gc"] - g2["contenido_gc"]), 2),
                "misma_hebra": g1["hebra"] == g2["hebra"],
            })

    # Genes unicos (no emparejados)
    all_prod_1 = set(idx_1.keys())
    all_prod_2 = set(idx_2.keys())
    prod_solo_1 = all_prod_1 - all_prod_2
    prod_solo_2 = all_prod_2 - all_prod_1

    unicos_1 = []
    for prod in prod_solo_1:
        for g in idx_1[prod]:
            unicos_1.append({"locus": g["locus_tag"], "nombre": g["nombre_gen"],
                             "producto": g["producto"], "longitud": g["longitud_pb"],
                             "inicio": g["inicio"], "fin": g["fin"],
                             "gc": g["contenido_gc"]})

    unicos_2 = []
    for prod in prod_solo_2:
        for g in idx_2[prod]:
            unicos_2.append({"locus": g["locus_tag"], "nombre": g["nombre_gen"],
                             "producto": g["producto"], "longitud": g["longitud_pb"],
                             "inicio": g["inicio"], "fin": g["fin"],
                             "gc": g["contenido_gc"]})

    # Hypothetical proteins
    hyp_1 = [g for g in genes_1 if "hypothetical" in g["producto"].lower()]
    hyp_2 = [g for g in genes_2 if "hypothetical" in g["producto"].lower()]

    total_1 = len(genes_1)
    total_2 = len(genes_2)

    # Ortologos con mayor diferencia de tamano
    ortologos_divergentes = sorted(ortologos, key=lambda x: x["dif_longitud"], reverse=True)[:20]

    # Estadisticas
    stats = {
        "total_ortologos": len(ortologos),
        "unicos_genoma_1": len(unicos_1),
        "unicos_genoma_2": len(unicos_2),
        "hipoteticos_genoma_1": len(hyp_1),
        "hipoteticos_genoma_2": len(hyp_2),
        "porcentaje_compartido_1": round(len(ortologos) / total_1 * 100, 1) if total_1 > 0 else 0,
        "porcentaje_compartido_2": round(len(ortologos) / total_2 * 100, 1) if total_2 > 0 else 0,
        "conservacion_hebra": sum(1 for o in ortologos if o["misma_hebra"]),
        "cambio_hebra": sum(1 for o in ortologos if not o["misma_hebra"]),
    }

    print(f"\n  Ortologos detectados: {stats['total_ortologos']:,}")
    print(f"  Unicos en {nombre_1}: {stats['unicos_genoma_1']:,}")
    print(f"  Unicos en {nombre_2}: {stats['unicos_genoma_2']:,}")
    print(f"  % compartido (G1): {stats['porcentaje_compartido_1']}%")
    print(f"  % compartido (G2): {stats['porcentaje_compartido_2']}%")

    return {
        "estadisticas": stats,
        "ortologos_muestra": [{"producto": o["producto"],
                                "gen_1_locus": o["gen_1"]["locus"], "gen_1_nombre": o["gen_1"]["nombre"],
                                "gen_2_locus": o["gen_2"]["locus"], "gen_2_nombre": o["gen_2"]["nombre"],
                                "dif_longitud": o["dif_longitud"], "dif_gc": o["dif_gc"],
                                "misma_hebra": o["misma_hebra"]}
                               for o in ortologos[:100]],
        "ortologos_divergentes": [{"producto": o["producto"],
                                    "gen_1_locus": o["gen_1"]["locus"],
                                    "gen_2_locus": o["gen_2"]["locus"],
                                    "longitud_1": o["gen_1"]["longitud"],
                                    "longitud_2": o["gen_2"]["longitud"],
                                    "dif_longitud": o["dif_longitud"],
                                    "dif_gc": o["dif_gc"]}
                                   for o in ortologos_divergentes],
        "unicos_genoma_1": unicos_1[:50],
        "unicos_genoma_2": unicos_2[:50],
        "total_unicos_1": len(unicos_1),
        "total_unicos_2": len(unicos_2),
        "_ortologos_full": ortologos  # para uso interno
    }


def analizar_mutaciones_ortologos(genes_1, genes_2, ortologos_data, nombre_1, nombre_2):
    """
    Analiza mutaciones en genes ortologos comparando secuencias de ADN.
    Para genes del mismo tamano: comparacion codon a codon.
    Clasifica mutaciones en sinonimas/no-sinonimas y transiciones/transversiones.
    """
    print("\n" + "=" * 70)
    print("ANALISIS DE MUTACIONES EN ORTOLOGOS")
    print("=" * 70)

    ortologos = ortologos_data.get("_ortologos_full", [])
    if not ortologos:
        print("  No hay ortologos para analizar")
        return {"error": "No hay ortologos"}

    # Crear indice de secuencias por locus_tag
    seq_idx_1 = {g["locus_tag"]: g["secuencia_adn"] for g in genes_1}
    seq_idx_2 = {g["locus_tag"]: g["secuencia_adn"] for g in genes_2}

    resultados = []
    total_snps = 0
    total_sinonimas = 0
    total_no_sinonimas = 0
    total_transiciones = 0
    total_transversiones = 0
    total_codones_comparados = 0
    identidades = []

    # Limitar a los primeros 500 ortologos para rendimiento
    for orto in ortologos[:500]:
        locus1 = orto["gen_1"]["locus"]
        locus2 = orto["gen_2"]["locus"]
        seq1 = seq_idx_1.get(locus1, "")
        seq2 = seq_idx_2.get(locus2, "")

        if not seq1 or not seq2:
            continue

        # Comparar las secuencias
        min_len = min(len(seq1), len(seq2))
        if min_len < 6:
            continue

        # Contar mismatches nucleotido a nucleotido
        matches = 0
        mismatches = 0
        for i in range(min_len):
            if seq1[i] == seq2[i]:
                matches += 1
            else:
                mismatches += 1

        # Agregar nucleotidos no cubiertos como diferencias
        len_diff = abs(len(seq1) - len(seq2))
        identity = round(matches / max(len(seq1), len(seq2)) * 100, 2)
        identidades.append(identity)

        # Analisis codon a codon (para genes del tamano similar)
        snps_gen = 0
        sinonimas = 0
        no_sinonimas = 0
        transiciones = 0
        transversiones = 0
        codones_comp = 0

        # Solo analizar codones si longitudes similares (diff < 10%)
        if len_diff <= max(len(seq1), len(seq2)) * 0.1:
            codon_len = (min_len // 3) * 3
            for i in range(0, codon_len, 3):
                c1 = seq1[i:i+3].upper()
                c2 = seq2[i:i+3].upper()
                if len(c1) != 3 or len(c2) != 3:
                    continue
                if not all(n in "ATGC" for n in c1 + c2):
                    continue

                codones_comp += 1
                if c1 != c2:
                    snps_gen += 1
                    # Traducir
                    aa1 = CODON_TABLE.get(c1, "X")
                    aa2 = CODON_TABLE.get(c2, "X")
                    if c1 in STOP_CODONS:
                        aa1 = "*"
                    if c2 in STOP_CODONS:
                        aa2 = "*"

                    if aa1 == aa2:
                        sinonimas += 1
                    else:
                        no_sinonimas += 1

                    # Tipo de sustitucion nucleotidica
                    for j in range(3):
                        n1 = c1[j]
                        n2 = c2[j]
                        if n1 != n2:
                            pair = frozenset([n1, n2])
                            if pair in [frozenset(["A", "G"]), frozenset(["C", "T"])]:
                                transiciones += 1
                            else:
                                transversiones += 1

        total_snps += snps_gen
        total_sinonimas += sinonimas
        total_no_sinonimas += no_sinonimas
        total_transiciones += transiciones
        total_transversiones += transversiones
        total_codones_comparados += codones_comp

        if mismatches > 0 or len_diff > 0:
            resultados.append({
                "producto": orto["producto"],
                "gen_1": locus1, "gen_2": locus2,
                "longitud_1": len(seq1), "longitud_2": len(seq2),
                "dif_longitud": len_diff,
                "snps_codones": snps_gen,
                "sinonimas": sinonimas, "no_sinonimas": no_sinonimas,
                "identidad_nt": identity,
                "codones_comparados": codones_comp,
            })

    # Ordenar por menor identidad (mas divergentes primero)
    resultados.sort(key=lambda x: x["identidad_nt"])

    # Calcular ANI basado en ortologos
    ani = round(statistics.mean(identidades), 2) if identidades else 0

    # Ka/Ks ratio (dN/dS)
    ka_ks = 0
    if total_sinonimas > 0 and total_no_sinonimas > 0:
        # Aproximacion simple
        ka_ks = round(total_no_sinonimas / total_sinonimas, 4)

    stats = {
        "ortologos_analizados": len(identidades),
        "ani_ortologos": ani,
        "total_snps_codones": total_snps,
        "mutaciones_sinonimas": total_sinonimas,
        "mutaciones_no_sinonimas": total_no_sinonimas,
        "ratio_ka_ks": ka_ks,
        "transiciones": total_transiciones,
        "transversiones": total_transversiones,
        "ratio_ti_tv": round(total_transiciones / total_transversiones, 2) if total_transversiones > 0 else 0,
        "codones_comparados": total_codones_comparados,
        "identidad_promedio": ani,
        "identidad_min": round(min(identidades), 2) if identidades else 0,
        "identidad_max": round(max(identidades), 2) if identidades else 0,
    }

    # Distribucion de identidades
    dist_identidad = {"100%": 0, "99-100%": 0, "95-99%": 0, "90-95%": 0, "80-90%": 0, "<80%": 0}
    for ident in identidades:
        if ident == 100: dist_identidad["100%"] += 1
        elif ident >= 99: dist_identidad["99-100%"] += 1
        elif ident >= 95: dist_identidad["95-99%"] += 1
        elif ident >= 90: dist_identidad["90-95%"] += 1
        elif ident >= 80: dist_identidad["80-90%"] += 1
        else: dist_identidad["<80%"] += 1

    print(f"\n  ANI (ortologos): {ani}%")
    print(f"  Total SNPs en codones: {total_snps:,}")
    print(f"  Sinonimas: {total_sinonimas:,} | No-sinonimas: {total_no_sinonimas:,}")
    print(f"  Ka/Ks: {ka_ks}")
    print(f"  Ti/Tv: {stats['ratio_ti_tv']}")

    return {
        "estadisticas": stats,
        "distribucion_identidad": dist_identidad,
        "genes_mas_divergentes": resultados[:30],
        "genes_identicos": [r for r in resultados if r["identidad_nt"] == 100][:20],
    }


def analizar_sintenia(ortologos_data, genes_1, genes_2, nombre_1, nombre_2):
    """
    Analiza conservacion del orden genico (sintenia).
    Detecta inversiones, translocaciones y bloques sintenicos.
    """
    print("\n" + "=" * 70)
    print("ANALISIS DE SINTENIA (ORDEN GENICO)")
    print("=" * 70)

    ortologos = ortologos_data.get("_ortologos_full", [])
    if not ortologos:
        return {"error": "No hay ortologos"}

    # Crear puntos para dot plot: (pos_genoma1, pos_genoma2)
    puntos = []
    for orto in ortologos:
        pos1 = (orto["gen_1"]["inicio"] + orto["gen_1"]["fin"]) / 2
        pos2 = (orto["gen_2"]["inicio"] + orto["gen_2"]["fin"]) / 2
        misma_hebra = orto["misma_hebra"]
        puntos.append({
            "x": round(pos1),
            "y": round(pos2),
            "misma_hebra": misma_hebra,
            "producto": orto["producto"][:50]
        })

    # Detectar bloques sintenicos (secuencias de genes en el mismo orden)
    # Ordenar ortologos por posicion en genoma 1
    orto_sorted = sorted(ortologos, key=lambda o: o["gen_1"]["inicio"])

    bloques = []
    bloque_actual = [orto_sorted[0]] if orto_sorted else []
    max_gap = 50000  # max 50kb gap entre genes del mismo bloque

    for i in range(1, len(orto_sorted)):
        prev = orto_sorted[i-1]
        curr = orto_sorted[i]

        # Verificar si estan en el mismo orden en ambos genomas
        prev_pos2 = prev["gen_2"]["inicio"]
        curr_pos2 = curr["gen_2"]["inicio"]
        same_direction = (curr_pos2 > prev_pos2) == (curr["gen_1"]["inicio"] > prev["gen_1"]["inicio"])

        gap1 = abs(curr["gen_1"]["inicio"] - prev["gen_1"]["fin"])
        gap2 = abs(curr["gen_2"]["inicio"] - prev["gen_2"]["fin"])

        if same_direction and gap1 < max_gap and gap2 < max_gap:
            bloque_actual.append(curr)
        else:
            if len(bloque_actual) >= 3:  # minimo 3 genes para un bloque
                bloques.append(bloque_actual)
            bloque_actual = [curr]

    if len(bloque_actual) >= 3:
        bloques.append(bloque_actual)

    # Procesar bloques
    bloques_info = []
    for bloque in bloques:
        inicio_1 = bloque[0]["gen_1"]["inicio"]
        fin_1 = bloque[-1]["gen_1"]["fin"]
        inicio_2 = min(b["gen_2"]["inicio"] for b in bloque)
        fin_2 = max(b["gen_2"]["fin"] for b in bloque)
        bloques_info.append({
            "genes": len(bloque),
            "inicio_g1": inicio_1, "fin_g1": fin_1,
            "inicio_g2": inicio_2, "fin_g2": fin_2,
            "longitud_g1": fin_1 - inicio_1,
            "longitud_g2": fin_2 - inicio_2,
        })

    # Detectar inversiones (genes en orden inverso en genoma 2)
    inversiones = 0
    for i in range(1, len(orto_sorted)):
        prev_pos2 = orto_sorted[i-1]["gen_2"]["inicio"]
        curr_pos2 = orto_sorted[i]["gen_2"]["inicio"]
        if curr_pos2 < prev_pos2:
            inversiones += 1

    # Detectar cambios de hebra
    cambios_hebra = sum(1 for o in ortologos if not o["misma_hebra"])

    total_genes_en_bloques = sum(b["genes"] for b in bloques_info)

    stats = {
        "total_puntos": len(puntos),
        "bloques_sintenicos": len(bloques_info),
        "genes_en_bloques": total_genes_en_bloques,
        "porcentaje_sintenico": round(total_genes_en_bloques / len(ortologos) * 100, 1) if ortologos else 0,
        "inversiones_detectadas": inversiones,
        "cambios_hebra": cambios_hebra,
    }

    print(f"\n  Bloques sintenicos: {stats['bloques_sintenicos']}")
    print(f"  Genes en bloques: {stats['genes_en_bloques']:,} ({stats['porcentaje_sintenico']}%)")
    print(f"  Inversiones: {stats['inversiones_detectadas']}")
    print(f"  Cambios de hebra: {stats['cambios_hebra']}")

    # Subsamplear puntos para el dot plot (max 2000 para rendimiento)
    if len(puntos) > 2000:
        puntos_plot = random.sample(puntos, 2000)
    else:
        puntos_plot = puntos

    return {
        "estadisticas": stats,
        "puntos_dotplot": puntos_plot,
        "bloques_sintenicos": bloques_info[:50],
    }


def detectar_islas_genomicas(genes_1, genes_2, ortologos_data, registro_1, registro_2, nombre_1, nombre_2):
    """
    Detecta islas genomicas: clusters de genes unicos + GC anormal.
    Las islas genomicas son regiones adquiridas por transferencia horizontal.
    """
    print("\n" + "=" * 70)
    print("DETECCION DE ISLAS GENOMICAS")
    print("=" * 70)

    gc_promedio_1 = gc_fraction(registro_1.seq) * 100
    gc_promedio_2 = gc_fraction(registro_2.seq) * 100

    # Marcar genes como compartidos o unicos
    ortologos = ortologos_data.get("_ortologos_full", [])
    locus_compartidos_1 = set(o["gen_1"]["locus"] for o in ortologos)
    locus_compartidos_2 = set(o["gen_2"]["locus"] for o in ortologos)

    def encontrar_islas(genes, locus_compartidos, gc_promedio, nombre):
        """Encuentra clusters de genes unicos consecutivos."""
        genes_sorted = sorted(genes, key=lambda g: g["inicio"])
        islas = []
        cluster = []

        for g in genes_sorted:
            es_unico = g["locus_tag"] not in locus_compartidos
            if es_unico:
                cluster.append(g)
            else:
                if len(cluster) >= 3:  # minimo 3 genes unicos consecutivos
                    islas.append(cluster)
                cluster = []

        if len(cluster) >= 3:
            islas.append(cluster)

        islas_info = []
        for isla in islas:
            inicio = isla[0]["inicio"]
            fin = isla[-1]["fin"]
            longitud = fin - inicio
            gc_isla = statistics.mean([g["contenido_gc"] for g in isla]) if isla else 0
            desviacion_gc = round(gc_isla - gc_promedio, 2)

            islas_info.append({
                "inicio": inicio,
                "fin": fin,
                "longitud_pb": longitud,
                "num_genes": len(isla),
                "gc_promedio": round(gc_isla, 2),
                "gc_genoma": round(gc_promedio, 2),
                "desviacion_gc": desviacion_gc,
                "posible_hgt": abs(desviacion_gc) > 3,  # >3% desviacion sugiere HGT
                "genes": [{"locus": g["locus_tag"], "nombre": g["nombre_gen"],
                           "producto": g["producto"][:60]} for g in isla[:10]],
            })

        return sorted(islas_info, key=lambda x: x["num_genes"], reverse=True)

    islas_1 = encontrar_islas(genes_1, locus_compartidos_1, gc_promedio_1, nombre_1)
    islas_2 = encontrar_islas(genes_2, locus_compartidos_2, gc_promedio_2, nombre_2)

    stats = {
        "islas_genoma_1": len(islas_1),
        "islas_genoma_2": len(islas_2),
        "genes_en_islas_1": sum(i["num_genes"] for i in islas_1),
        "genes_en_islas_2": sum(i["num_genes"] for i in islas_2),
        "pb_en_islas_1": sum(i["longitud_pb"] for i in islas_1),
        "pb_en_islas_2": sum(i["longitud_pb"] for i in islas_2),
        "posibles_hgt_1": sum(1 for i in islas_1 if i["posible_hgt"]),
        "posibles_hgt_2": sum(1 for i in islas_2 if i["posible_hgt"]),
    }

    print(f"\n  Islas en {nombre_1}: {stats['islas_genoma_1']} ({stats['genes_en_islas_1']} genes)")
    print(f"  Islas en {nombre_2}: {stats['islas_genoma_2']} ({stats['genes_en_islas_2']} genes)")
    print(f"  Posible HGT en {nombre_1}: {stats['posibles_hgt_1']}")
    print(f"  Posible HGT en {nombre_2}: {stats['posibles_hgt_2']}")

    return {
        "estadisticas": stats,
        "islas_genoma_1": islas_1[:20],
        "islas_genoma_2": islas_2[:20],
    }


def comparar_categorias_funcionales(genes_1, genes_2, nombre_1, nombre_2):
    """
    Compara el perfil funcional de ambos genomas basado en anotaciones.
    """
    print("\n" + "=" * 70)
    print("COMPARACION DE PERFIL FUNCIONAL")
    print("=" * 70)

    def categorizar(genes):
        conteo = Counter()
        for g in genes:
            prod = g["producto"].lower()
            categorizada = False
            for cat, keywords in CATEGORIAS_FUNCIONALES.items():
                if any(kw.lower() in prod for kw in keywords):
                    conteo[cat] += 1
                    categorizada = True
                    break
            if not categorizada:
                conteo["Otros"] += 1
        return dict(conteo)

    cats_1 = categorizar(genes_1)
    cats_2 = categorizar(genes_2)

    # Todas las categorias
    todas = sorted(set(list(cats_1.keys()) + list(cats_2.keys())))

    comparacion = []
    for cat in todas:
        c1 = cats_1.get(cat, 0)
        c2 = cats_2.get(cat, 0)
        comparacion.append({
            "categoria": cat,
            "genoma_1": c1,
            "genoma_2": c2,
            "diferencia": c2 - c1,
            "porcentaje_1": round(c1 / len(genes_1) * 100, 1) if genes_1 else 0,
            "porcentaje_2": round(c2 / len(genes_2) * 100, 1) if genes_2 else 0,
        })

    # Ordenar por mayor diferencia absoluta
    comparacion.sort(key=lambda x: abs(x["diferencia"]), reverse=True)

    for c in comparacion[:10]:
        print(f"  {c['categoria']:<35} G1: {c['genoma_1']:>5} ({c['porcentaje_1']}%)  "
              f"G2: {c['genoma_2']:>5} ({c['porcentaje_2']}%)  Dif: {c['diferencia']:>+4}")

    return {
        "categorias": comparacion,
        "total_categorizado_1": sum(cats_1.values()),
        "total_categorizado_2": sum(cats_2.values()),
    }


def comparar_trna_rrna(registro_1, registro_2, nombre_1, nombre_2):
    """Compara tRNA y rRNA entre los genomas."""
    print("\n" + "=" * 70)
    print("COMPARACION DE tRNA Y rRNA")
    print("=" * 70)

    feat_1 = extraer_features_no_cds(registro_1)
    feat_2 = extraer_features_no_cds(registro_2)

    # tRNA por aminoacido
    trna_aa_1 = Counter(t["aminoacido"] for t in feat_1["tRNA"] if t.get("aminoacido"))
    trna_aa_2 = Counter(t["aminoacido"] for t in feat_2["tRNA"] if t.get("aminoacido"))

    todos_aa = sorted(set(list(trna_aa_1.keys()) + list(trna_aa_2.keys())))
    trna_comparacion = []
    for aa in todos_aa:
        c1 = trna_aa_1.get(aa, 0)
        c2 = trna_aa_2.get(aa, 0)
        trna_comparacion.append({"aminoacido": aa, "genoma_1": c1, "genoma_2": c2, "diferencia": c2 - c1})

    # rRNA por subtipo
    rrna_sub_1 = Counter(r["subtipo"] for r in feat_1["rRNA"] if r.get("subtipo"))
    rrna_sub_2 = Counter(r["subtipo"] for r in feat_2["rRNA"] if r.get("subtipo"))

    todos_sub = sorted(set(list(rrna_sub_1.keys()) + list(rrna_sub_2.keys())))
    rrna_comparacion = []
    for sub in todos_sub:
        c1 = rrna_sub_1.get(sub, 0)
        c2 = rrna_sub_2.get(sub, 0)
        rrna_comparacion.append({"subtipo": sub, "genoma_1": c1, "genoma_2": c2, "diferencia": c2 - c1})

    stats = {
        "total_trna_1": len(feat_1["tRNA"]),
        "total_trna_2": len(feat_2["tRNA"]),
        "total_rrna_1": len(feat_1["rRNA"]),
        "total_rrna_2": len(feat_2["rRNA"]),
        "total_ncrna_1": len(feat_1.get("ncRNA", [])),
        "total_ncrna_2": len(feat_2.get("ncRNA", [])),
        "tipos_trna_1": len(trna_aa_1),
        "tipos_trna_2": len(trna_aa_2),
    }

    print(f"\n  tRNA: {nombre_1}={stats['total_trna_1']}, {nombre_2}={stats['total_trna_2']}")
    print(f"  rRNA: {nombre_1}={stats['total_rrna_1']}, {nombre_2}={stats['total_rrna_2']}")
    print(f"  ncRNA: {nombre_1}={stats['total_ncrna_1']}, {nombre_2}={stats['total_ncrna_2']}")

    return {
        "estadisticas": stats,
        "trna_por_aminoacido": trna_comparacion,
        "rrna_por_subtipo": rrna_comparacion,
    }


def analizar_resistencia_antibiotica(genes_1, genes_2, nombre_1, nombre_2):
    """Detecta y compara genes de resistencia antibiotica."""
    print("\n" + "=" * 70)
    print("ANALISIS DE RESISTENCIA ANTIBIOTICA")
    print("=" * 70)

    def detectar_resistencia(genes):
        resultado = {}
        for categoria, keywords in GENES_RESISTENCIA.items():
            genes_cat = []
            for g in genes:
                prod = g["producto"].lower()
                if any(kw.lower() in prod for kw in keywords):
                    genes_cat.append({
                        "locus": g["locus_tag"],
                        "nombre": g["nombre_gen"],
                        "producto": g["producto"],
                    })
            resultado[categoria] = {
                "total": len(genes_cat),
                "genes": genes_cat[:10]
            }
        return resultado

    res_1 = detectar_resistencia(genes_1)
    res_2 = detectar_resistencia(genes_2)

    # Comparacion
    comparacion = []
    for cat in GENES_RESISTENCIA.keys():
        c1 = res_1[cat]["total"]
        c2 = res_2[cat]["total"]
        comparacion.append({
            "categoria": cat,
            "genoma_1": c1,
            "genoma_2": c2,
            "diferencia": c2 - c1,
            "genes_1": res_1[cat]["genes"][:5],
            "genes_2": res_2[cat]["genes"][:5],
        })

    total_1 = sum(r["total"] for r in res_1.values())
    total_2 = sum(r["total"] for r in res_2.values())

    stats = {
        "total_genes_resistencia_1": total_1,
        "total_genes_resistencia_2": total_2,
        "categorias_con_genes_1": sum(1 for r in res_1.values() if r["total"] > 0),
        "categorias_con_genes_2": sum(1 for r in res_2.values() if r["total"] > 0),
    }

    print(f"\n  Total genes resistencia {nombre_1}: {total_1}")
    print(f"  Total genes resistencia {nombre_2}: {total_2}")
    for c in comparacion:
        if c["genoma_1"] > 0 or c["genoma_2"] > 0:
            print(f"    {c['categoria']:<25} G1: {c['genoma_1']:>3}  G2: {c['genoma_2']:>3}")

    return {
        "estadisticas": stats,
        "comparacion": comparacion,
    }


def calcular_gc_ventana(registro_1, registro_2, nombre_1, nombre_2, ventana=50000):
    """
    Calcula GC% en ventanas deslizantes para ambos genomas.
    Permite detectar regiones con GC anormal (posible HGT).
    """
    print("\n" + "=" * 70)
    print("ANALISIS DE GC REGIONAL (VENTANA DESLIZANTE)")
    print("=" * 70)

    def gc_sliding(registro, step=None):
        seq = str(registro.seq).upper()
        if step is None:
            step = ventana // 2
        puntos = []
        for i in range(0, len(seq) - ventana, step):
            fragmento = seq[i:i + ventana]
            gc = (fragmento.count("G") + fragmento.count("C")) / len(fragmento) * 100
            puntos.append({
                "posicion": i + ventana // 2,
                "gc": round(gc, 2)
            })
        return puntos

    gc_1 = gc_sliding(registro_1)
    gc_2 = gc_sliding(registro_2)

    # Detectar regiones con GC anormal (>2 std dev de la media)
    gc_vals_1 = [p["gc"] for p in gc_1]
    gc_vals_2 = [p["gc"] for p in gc_2]

    mean_1 = statistics.mean(gc_vals_1) if gc_vals_1 else 0
    std_1 = statistics.stdev(gc_vals_1) if len(gc_vals_1) > 1 else 0
    mean_2 = statistics.mean(gc_vals_2) if gc_vals_2 else 0
    std_2 = statistics.stdev(gc_vals_2) if len(gc_vals_2) > 1 else 0

    anomalas_1 = [p for p in gc_1 if abs(p["gc"] - mean_1) > 2 * std_1] if std_1 > 0 else []
    anomalas_2 = [p for p in gc_2 if abs(p["gc"] - mean_2) > 2 * std_2] if std_2 > 0 else []

    # Subsamplear para el frontend (max 200 puntos)
    max_pts = 200
    if len(gc_1) > max_pts:
        step = len(gc_1) // max_pts
        gc_1_plot = gc_1[::step][:max_pts]
    else:
        gc_1_plot = gc_1

    if len(gc_2) > max_pts:
        step = len(gc_2) // max_pts
        gc_2_plot = gc_2[::step][:max_pts]
    else:
        gc_2_plot = gc_2

    stats = {
        "ventana_pb": ventana,
        "media_gc_1": round(mean_1, 2),
        "media_gc_2": round(mean_2, 2),
        "std_gc_1": round(std_1, 2),
        "std_gc_2": round(std_2, 2),
        "regiones_anomalas_1": len(anomalas_1),
        "regiones_anomalas_2": len(anomalas_2),
    }

    print(f"\n  Ventana: {ventana:,} pb")
    print(f"  GC media {nombre_1}: {stats['media_gc_1']}% (std: {stats['std_gc_1']}%)")
    print(f"  GC media {nombre_2}: {stats['media_gc_2']}% (std: {stats['std_gc_2']}%)")
    print(f"  Regiones anomalas {nombre_1}: {stats['regiones_anomalas_1']}")
    print(f"  Regiones anomalas {nombre_2}: {stats['regiones_anomalas_2']}")

    return {
        "estadisticas": stats,
        "gc_genoma_1": gc_1_plot,
        "gc_genoma_2": gc_2_plot,
        "anomalas_genoma_1": anomalas_1[:20],
        "anomalas_genoma_2": anomalas_2[:20],
    }


def estimar_ani_fragmentos(registro_1, registro_2, nombre_1, nombre_2, n_fragmentos=200, tam_fragmento=500):
    """
    Estima ANI (Average Nucleotide Identity) usando fragmentos aleatorios.
    Busca matches exactos de fragmentos en el otro genoma.
    """
    print("\n" + "=" * 70)
    print("ESTIMACION DE ANI (IDENTIDAD NUCLEOTIDICA PROMEDIO)")
    print("=" * 70)

    seq1 = str(registro_1.seq).upper()
    seq2 = str(registro_2.seq).upper()

    identidades = []
    matches_found = 0

    random.seed(42)  # reproducibilidad
    for _ in range(n_fragmentos):
        # Fragmento aleatorio de genoma 1
        pos = random.randint(0, len(seq1) - tam_fragmento)
        frag = seq1[pos:pos + tam_fragmento]

        # Buscar seed (20-mer) en genoma 2
        seed = frag[:20]
        idx = seq2.find(seed)
        if idx >= 0:
            # Comparar fragmento completo
            frag2 = seq2[idx:idx + tam_fragmento]
            if len(frag2) == tam_fragmento:
                matches = sum(1 for a, b in zip(frag, frag2) if a == b)
                ident = matches / tam_fragmento * 100
                identidades.append(ident)
                matches_found += 1

    ani = round(statistics.mean(identidades), 2) if identidades else 0
    cobertura = round(matches_found / n_fragmentos * 100, 1)

    print(f"\n  Fragmentos analizados: {n_fragmentos}")
    print(f"  Matches encontrados: {matches_found} ({cobertura}%)")
    print(f"  ANI estimado: {ani}%")

    # Clasificacion taxonomica basada en ANI
    if ani >= 96:
        clasificacion = "Misma especie (ANI >= 96%)"
    elif ani >= 90:
        clasificacion = "Mismo genero (ANI 90-96%)"
    elif ani >= 80:
        clasificacion = "Misma familia (ANI 80-90%)"
    else:
        clasificacion = "Familias diferentes (ANI < 80%)"

    return {
        "ani_estimado": ani,
        "fragmentos_total": n_fragmentos,
        "fragmentos_match": matches_found,
        "cobertura": cobertura,
        "clasificacion": clasificacion,
        "tamano_fragmento": tam_fragmento,
    }


# =============================================================================
# INTERPRETACION IA
# =============================================================================

def obtener_interpretacion_ia(resultados_completos, nombre_1, nombre_2):
    """Obtiene interpretacion biologica via Gemini."""
    print("\n[INFO] Consultando IA para interpretacion biologica...")
    try:
        sys.path.insert(0, DIRECTORIO_PROYECTO)
        from backend.gemini_client import consultar_gemini

        metricas = resultados_completos.get("metricas_generales", {})
        ec = metricas.get("ecoli", {})
        sal = metricas.get("salmonella", {})
        vir = resultados_completos.get("genes_virulencia", {})
        orto = resultados_completos.get("ortologos", {}).get("estadisticas", {})
        mut = resultados_completos.get("mutaciones", {}).get("estadisticas", {})
        ani_data = resultados_completos.get("ani", {})

        resumen = f"""Comparacion genomica avanzada entre {nombre_1} y {nombre_2}:
METRICAS BASICAS:
- {nombre_1}: {ec.get('longitud_genoma_pb', 0):,} pb, {ec.get('total_genes_cds', 0)} CDS, GC {ec.get('contenido_gc_porcentaje', 0)}%
- {nombre_2}: {sal.get('longitud_genoma_pb', 0):,} pb, {sal.get('total_genes_cds', 0)} CDS, GC {sal.get('contenido_gc_porcentaje', 0)}%
ORTOLOGOS:
- Genes compartidos: {orto.get('total_ortologos', 0)}, Unicos G1: {orto.get('unicos_genoma_1', 0)}, Unicos G2: {orto.get('unicos_genoma_2', 0)}
MUTACIONES:
- ANI (ortologos): {mut.get('ani_ortologos', 0)}%, Ka/Ks: {mut.get('ratio_ka_ks', 0)}, Ti/Tv: {mut.get('ratio_ti_tv', 0)}
- SNPs: {mut.get('total_snps_codones', 0)}, Sinonimas: {mut.get('mutaciones_sinonimas', 0)}, No-sinonimas: {mut.get('mutaciones_no_sinonimas', 0)}
ANI (fragmentos): {ani_data.get('ani_estimado', 0)}% - {ani_data.get('clasificacion', '')}
VIRULENCIA: {nombre_1}={vir.get('ecoli', {}).get('total', 0)}, {nombre_2}={vir.get('salmonella', {}).get('total', 0)}"""

        prompt = f"""Eres un experto en genomica comparativa y evolucion bacteriana. Analiza esta comparacion avanzada y da una interpretacion biologica profunda en espanol (5-7 parrafos).

{resumen}

Incluye:
1. Relacion evolutiva entre ambos organismos (basado en ANI y ortologos)
2. Significado del ratio Ka/Ks (seleccion purificadora, neutral o positiva)
3. Interpretacion de las mutaciones (transiciones vs transversiones)
4. Diferencias en capacidad patogenica
5. Hipotesis sobre transferencia horizontal de genes
6. Implicaciones para la cepa cultivada/derivada (si aplica)
7. Conclusiones principales

Responde de forma directa y cientifica sin introducciones."""

        respuesta, error = consultar_gemini(prompt)
        if respuesta:
            print("[INFO] Interpretacion IA obtenida")
            return respuesta
        else:
            print(f"[INFO] IA no disponible: {error}")
            return None
    except Exception as e:
        print(f"[INFO] Error IA: {e}")
        return None


def generar_resumen_comparativo(resultados, nombre_1, nombre_2):
    """Genera resumen con todas las diferencias clave."""
    metricas = resultados.get("metricas_generales", {})
    orto = resultados.get("ortologos", {}).get("estadisticas", {})
    mut = resultados.get("mutaciones", {}).get("estadisticas", {})
    islas = resultados.get("islas_genomicas", {}).get("estadisticas", {})

    resumen = {
        "titulo": f"Comparacion: {nombre_1} vs {nombre_2}",
        "diferencias_clave": [
            {"aspecto": "Tamano del genoma",
             "observacion": f"Diferencia de {abs(metricas.get('diferencias', {}).get('longitud_pb', 0)):,} pb",
             "significado": "ADN adicional puede contener genes de adaptacion"},
            {"aspecto": "Genes compartidos",
             "observacion": f"{orto.get('total_ortologos', 0):,} ortologos, {orto.get('unicos_genoma_1', 0):,} unicos G1, {orto.get('unicos_genoma_2', 0):,} unicos G2",
             "significado": "Genes unicos definen las capacidades especificas de cada cepa"},
            {"aspecto": "Mutaciones en ortologos",
             "observacion": f"ANI: {mut.get('ani_ortologos', 0)}%, Ka/Ks: {mut.get('ratio_ka_ks', 0)}",
             "significado": "Ka/Ks < 1 = seleccion purificadora, > 1 = seleccion positiva"},
            {"aspecto": "Islas genomicas",
             "observacion": f"G1: {islas.get('islas_genoma_1', 0)} islas, G2: {islas.get('islas_genoma_2', 0)} islas",
             "significado": "Regiones probablemente adquiridas por transferencia horizontal"},
        ],
        "conclusiones": [
            "La comparacion revela diferencias en contenido genico y organizacion genomica",
            "Los genes unicos definen las capacidades adaptativas especificas de cada organismo",
            "El analisis de mutaciones indica el grado de divergencia evolutiva",
        ]
    }
    return resumen


# =============================================================================
# EXPORTACION
# =============================================================================

def exportar_comparacion_json(datos, nombre_archivo):
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.json")
    # Limpiar datos internos antes de exportar
    datos_export = json.loads(json.dumps(datos, default=str))
    if "ortologos" in datos_export and "_ortologos_full" in datos_export["ortologos"]:
        del datos_export["ortologos"]["_ortologos_full"]
    with open(ruta, 'w', encoding='utf-8') as f:
        json.dump(datos_export, f, indent=2, ensure_ascii=False)
    print(f"       [OK] {nombre_archivo}.json")


def exportar_comparacion_csv(metricas, nombre_archivo, nombre_1, nombre_2):
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
        ]
        for fila in filas:
            writer.writerow(fila)
    print(f"       [OK] {nombre_archivo}.csv")


def exportar_genes_virulencia_csv(virulencia_datos, nombre_archivo):
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")
    with open(ruta, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Locus_tag', 'Nombre_gen', 'Producto'])
        for gen in virulencia_datos["ejemplos"]:
            writer.writerow([gen['locus'], gen['gen'], gen['producto']])
    print(f"       [OK] {nombre_archivo}.csv")


# =============================================================================
# FUNCION PRINCIPAL
# =============================================================================

def main():
    print("\n" + "=" * 70)
    print("COMPARACION AVANZADA DE GENOMAS BACTERIANOS")
    print("=" * 70)
    print(f"\n  Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Genoma 1: {BASENAME_1}")
    print(f"  Genoma 2: {BASENAME_2}")

    crear_directorios()

    archivos_ok, faltantes = verificar_archivos()
    if not archivos_ok:
        print("[ERROR] Faltan archivos GenBank:")
        for f in faltantes:
            print(f"        - {f}")
        return

    # PASO 1: Cargar genomas
    print("\n[PASO 1/8] Cargando archivos GenBank...")
    registro_1 = cargar_genbank(ARCHIVO_GB_1, BASENAME_1)
    registro_2 = cargar_genbank(ARCHIVO_GB_2, BASENAME_2)
    nombre_1 = obtener_nombre_organismo(registro_1, BASENAME_1)
    nombre_2 = obtener_nombre_organismo(registro_2, BASENAME_2)
    print(f"       Organismo 1: {nombre_1}")
    print(f"       Organismo 2: {nombre_2}")

    # PASO 2: Extraer genes
    print("\n[PASO 2/8] Extrayendo genes...")
    genes_1 = extraer_genes_con_info(registro_1, nombre_1)
    genes_2 = extraer_genes_con_info(registro_2, nombre_2)
    print(f"       {nombre_1}: {len(genes_1):,} CDS")
    print(f"       {nombre_2}: {len(genes_2):,} CDS")

    # PASO 3: Comparaciones basicas
    print("\n[PASO 3/8] Comparaciones basicas...")
    metricas = comparar_metricas_generales(registro_1, registro_2, genes_1, genes_2, nombre_1, nombre_2)
    virulencia = comparar_virulencia(genes_1, genes_2, nombre_1, nombre_2)
    distancias = comparar_distancias_intergenicas(genes_1, genes_2, nombre_1, nombre_2)
    uso_codones = comparar_uso_codones(registro_1, registro_2, genes_1, genes_2, nombre_1, nombre_2)
    dist_tamanos = comparar_distribucion_tamanos(genes_1, genes_2, nombre_1, nombre_2)
    dist_gc = comparar_distribucion_gc(genes_1, genes_2, nombre_1, nombre_2)

    # PASO 4: Deteccion de ortologos
    print("\n[PASO 4/8] Detectando ortologos...")
    ortologos = detectar_ortologos(genes_1, genes_2, nombre_1, nombre_2)

    # PASO 5: Analisis de mutaciones
    print("\n[PASO 5/8] Analizando mutaciones en ortologos...")
    mutaciones = analizar_mutaciones_ortologos(genes_1, genes_2, ortologos, nombre_1, nombre_2)

    # PASO 6: Analisis avanzados
    print("\n[PASO 6/8] Analisis avanzados...")
    sintenia = analizar_sintenia(ortologos, genes_1, genes_2, nombre_1, nombre_2)
    islas = detectar_islas_genomicas(genes_1, genes_2, ortologos, registro_1, registro_2, nombre_1, nombre_2)
    funcional = comparar_categorias_funcionales(genes_1, genes_2, nombre_1, nombre_2)
    trna_rrna = comparar_trna_rrna(registro_1, registro_2, nombre_1, nombre_2)
    resistencia = analizar_resistencia_antibiotica(genes_1, genes_2, nombre_1, nombre_2)

    # PASO 7: GC regional y ANI
    print("\n[PASO 7/8] GC regional y ANI...")
    gc_regional = calcular_gc_ventana(registro_1, registro_2, nombre_1, nombre_2)
    ani = estimar_ani_fragmentos(registro_1, registro_2, nombre_1, nombre_2)

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
        "ortologos": ortologos,
        "mutaciones": mutaciones,
        "sintenia": sintenia,
        "islas_genomicas": islas,
        "categorias_funcionales": funcional,
        "trna_rrna": trna_rrna,
        "resistencia_antibiotica": resistencia,
        "gc_regional": gc_regional,
        "ani": ani,
    }

    # Resumen
    resumen = generar_resumen_comparativo(resultados_completos, nombre_1, nombre_2)
    resultados_completos["resumen_interpretativo"] = resumen

    # Interpretacion IA
    interpretacion = obtener_interpretacion_ia(resultados_completos, nombre_1, nombre_2)
    if interpretacion:
        resultados_completos["interpretacion_ia"] = interpretacion

    # PASO 8: Exportar
    print("\n[PASO 8/8] Exportando resultados...")
    nombre_comparacion = f"comparacion_{BASENAME_1}_vs_{BASENAME_2}"
    exportar_comparacion_json(resultados_completos, nombre_comparacion)
    exportar_comparacion_csv(metricas, f"metricas_{nombre_comparacion}", nombre_1, nombre_2)
    exportar_genes_virulencia_csv(virulencia["salmonella"], f"genes_virulencia_{BASENAME_2}")

    print("\n" + "=" * 70)
    print("[OK] COMPARACION AVANZADA COMPLETADA")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()

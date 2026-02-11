#!/usr/bin/env python3
"""
GenomeHub - Analisis Completo de Proteinas
==========================================
Analiza el proteoma completo de un genoma bacteriano:
- Estructura Primaria: composicion, propiedades fisicoquimicas, funciones, mutaciones
- Estructura Secundaria: prediccion helix/sheet/turn (Chou-Fasman via BioPython)
- Estructura Terciaria: consulta PDB/AlphaFold, analisis de cisteinas
- Estructura Cuaternaria: deteccion de complejos multi-subunidad

Uso: python analisis_proteinas.py <genome_basename>
"""

import os
import sys
import json
import csv
import gc
import re
import statistics
from datetime import datetime
from collections import Counter, defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False
    print("[AVISO] Modulo 'requests' no instalado. Consultas a PDB/AlphaFold deshabilitadas.")


# =============================================================================
# CONFIGURACION
# =============================================================================

DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # backend/
DIRECTORIO_PROYECTO = os.path.dirname(DIRECTORIO_BASE)  # raiz del proyecto
GENOME_BASENAME = sys.argv[1] if len(sys.argv) >= 2 else None

if not GENOME_BASENAME:
    print("[ERROR] Uso: python analisis_proteinas.py <genome_basename>")
    sys.exit(1)

# Buscar archivo GenBank
RUTA_DATOS = os.path.join(DIRECTORIO_PROYECTO, "datos", "crudo")
archivo_gb = None
for ext in [".gb", ".gbk", ".genbank"]:
    ruta = os.path.join(RUTA_DATOS, GENOME_BASENAME + ext)
    if os.path.exists(ruta):
        archivo_gb = ruta
        break

if not archivo_gb:
    for f in os.listdir(RUTA_DATOS):
        if GENOME_BASENAME in f and f.endswith((".gb", ".gbk", ".genbank")):
            archivo_gb = os.path.join(RUTA_DATOS, f)
            break

if not archivo_gb:
    print(f"[ERROR] No se encontro archivo GenBank para '{GENOME_BASENAME}' en {RUTA_DATOS}")
    sys.exit(1)

# Directorios de salida
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados")
RUTA_TABLAS = os.path.join(RUTA_RESULTADOS, GENOME_BASENAME, "tablas")
RUTA_FIGURAS = os.path.join(RUTA_RESULTADOS, GENOME_BASENAME, "figuras")

for d in [RUTA_TABLAS, RUTA_FIGURAS]:
    os.makedirs(d, exist_ok=True)

# Cache de consultas PDB/AlphaFold
RUTA_CACHE = os.path.join(RUTA_RESULTADOS, GENOME_BASENAME, "cache_proteinas.json")

# Estilo de graficos
try:
    plt.style.use('seaborn-v0_8-whitegrid')
except OSError:
    try:
        plt.style.use('seaborn-whitegrid')
    except OSError:
        plt.style.use('ggplot')

COLORES = {
    'primario': '#10B981',
    'secundario': '#A23B72',
    'terciario': '#F18F01',
    'cuaternario': '#C73E1D',
    'azul': '#3B82F6',
    'cyan': '#06B6D4',
    'gris': '#6C757D',
    'verde': '#3A7D44'
}

# Aminoacidos hidrofobicos para deteccion de peptido senal
AA_HIDROFOBICOS = set('AVILMFWP')

# Aminoacidos y sus propiedades
AA_NOMBRES = {
    'A': 'Alanina', 'R': 'Arginina', 'N': 'Asparagina', 'D': 'Acido aspartico',
    'C': 'Cisteina', 'E': 'Acido glutamico', 'Q': 'Glutamina', 'G': 'Glicina',
    'H': 'Histidina', 'I': 'Isoleucina', 'L': 'Leucina', 'K': 'Lisina',
    'M': 'Metionina', 'F': 'Fenilalanina', 'P': 'Prolina', 'S': 'Serina',
    'T': 'Treonina', 'W': 'Triptofano', 'Y': 'Tirosina', 'V': 'Valina'
}

# Genes criticos para analisis de mutaciones (E. coli y bacterias patogenas)
GENES_CRITICOS = {
    "gyrA": {
        "descripcion": "DNA girasa subunidad A - diana de fluoroquinolonas",
        "posiciones_clave": {83: "Ser->Leu/Trp (resistencia a ciprofloxacina)", 87: "Asp->Asn/Gly (resistencia a fluoroquinolonas)"},
        "longitud_esperada_aa": 875
    },
    "gyrB": {
        "descripcion": "DNA girasa subunidad B - diana de novobiocina",
        "posiciones_clave": {426: "Asp->Asn (resistencia)", 447: "Lys->Glu (resistencia)"},
        "longitud_esperada_aa": 804
    },
    "rpoB": {
        "descripcion": "RNA polimerasa subunidad beta - diana de rifampicina",
        "posiciones_clave": {516: "Asp->Val (resistencia a rifampicina)", 526: "His->Tyr/Asp (resistencia)", 531: "Ser->Leu (resistencia mas comun)"},
        "longitud_esperada_aa": 1342
    },
    "recA": {
        "descripcion": "Recombinasa - reparacion de ADN, respuesta SOS",
        "posiciones_clave": {72: "Glu clave para actividad ATPasa"},
        "longitud_esperada_aa": 353
    },
    "parC": {
        "descripcion": "Topoisomerasa IV subunidad A - diana de fluoroquinolonas",
        "posiciones_clave": {80: "Ser->Ile/Arg (resistencia a fluoroquinolonas)", 84: "Glu->Val/Gly (resistencia)"},
        "longitud_esperada_aa": 752
    },
    "parE": {
        "descripcion": "Topoisomerasa IV subunidad B",
        "posiciones_clave": {458: "Ser->Ala (resistencia)"},
        "longitud_esperada_aa": 630
    },
    "folA": {
        "descripcion": "Dihidrofolato reductasa - diana de trimetoprim",
        "posiciones_clave": {28: "Ile->Val (resistencia a trimetoprim)", 94: "Leu->Pro (resistencia)"},
        "longitud_esperada_aa": 159
    },
    "ampC": {
        "descripcion": "Beta-lactamasa cromosomal - resistencia a beta-lactamicos",
        "posiciones_clave": {},
        "longitud_esperada_aa": 381
    }
}


# =============================================================================
# FUNCIONES DE UTILIDAD
# =============================================================================

def fmt_num(n, decimals=0):
    """Formatea numeros con separador de miles."""
    if isinstance(n, float):
        return f"{n:,.{decimals}f}"
    return f"{n:,}"


def guardar_figura(fig, nombre):
    """Guarda una figura PNG en el directorio de figuras."""
    ruta = os.path.join(RUTA_FIGURAS, f"{nombre}.png")
    fig.savefig(ruta, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"  [OK] Figura: {nombre}.png")
    plt.close(fig)
    plt.close('all')
    gc.collect()


def cargar_cache():
    """Carga cache de consultas previas."""
    if os.path.exists(RUTA_CACHE):
        with open(RUTA_CACHE, 'r', encoding='utf-8') as f:
            return json.load(f)
    return {"pdb": {}, "alphafold": {}}


def guardar_cache(cache):
    """Guarda cache de consultas."""
    with open(RUTA_CACHE, 'w', encoding='utf-8') as f:
        json.dump(cache, f, ensure_ascii=False, indent=2)


def limpiar_secuencia_proteina(seq_str):
    """Limpia secuencia de proteina removiendo * (stop) y caracteres invalidos."""
    seq_str = seq_str.replace('*', '').replace('X', '').replace('J', '').replace('B', '').replace('Z', '')
    return seq_str.strip()


# =============================================================================
# 1. ESTRUCTURA PRIMARIA
# =============================================================================

def extraer_proteinas(registro):
    """
    Extrae y traduce todas las proteinas CDS del genoma.
    """
    print("\n" + "=" * 70)
    print("EXTRACCION Y TRADUCCION DE PROTEINAS")
    print("=" * 70)

    secuencia_genoma = registro.seq
    proteinas = []
    errores = 0

    for feature in registro.features:
        if feature.type != "CDS":
            continue

        try:
            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
            nombre_gen = feature.qualifiers.get("gene", [""])[0]
            producto = feature.qualifiers.get("product", [""])[0]
            proteina_id = feature.qualifiers.get("protein_id", [""])[0]

            inicio = int(feature.location.start)
            fin = int(feature.location.end)
            hebra = "+" if feature.location.strand == 1 else "-"
            longitud_pb = fin - inicio

            # Traducir secuencia de DNA a proteina
            secuencia_dna = feature.location.extract(secuencia_genoma)

            # Intentar obtener la traduccion anotada primero
            translation = feature.qualifiers.get("translation", [None])[0]
            if translation:
                secuencia_proteina = str(translation)
            else:
                # Traducir manualmente
                try:
                    secuencia_proteina = str(secuencia_dna.translate(table=11, to_stop=True))
                except Exception:
                    secuencia_proteina = str(secuencia_dna.translate(to_stop=True))

            secuencia_limpia = limpiar_secuencia_proteina(secuencia_proteina)
            if len(secuencia_limpia) < 5:
                continue

            proteinas.append({
                "locus_tag": locus_tag,
                "nombre_gen": nombre_gen,
                "producto": producto,
                "proteina_id": proteina_id,
                "inicio": inicio,
                "fin": fin,
                "hebra": hebra,
                "longitud_pb": longitud_pb,
                "longitud_aa": len(secuencia_limpia),
                "secuencia_proteina": secuencia_limpia
            })

        except Exception as e:
            errores += 1
            if errores <= 5:
                print(f"  [AVISO] Error traduciendo CDS: {e}")

    print(f"[OK] {len(proteinas):,} proteinas extraidas y traducidas")
    if errores > 0:
        print(f"  [AVISO] {errores} CDS no pudieron ser traducidos")

    return proteinas


def analizar_composicion_aminoacidos(proteinas):
    """
    Calcula la composicion de aminoacidos del proteoma completo.
    """
    print("\n" + "=" * 70)
    print("COMPOSICION DE AMINOACIDOS DEL PROTEOMA")
    print("=" * 70)

    aminoacidos_orden = list('ACDEFGHIKLMNPQRSTVWY')
    conteo_global = Counter()
    total_aa = 0

    for p in proteinas:
        for aa in p["secuencia_proteina"]:
            if aa in aminoacidos_orden:
                conteo_global[aa] += 1
                total_aa += 1

    composicion = {}
    for aa in aminoacidos_orden:
        porcentaje = (conteo_global[aa] / total_aa * 100) if total_aa > 0 else 0
        composicion[aa] = round(porcentaje, 2)
        print(f"  {aa} ({AA_NOMBRES.get(aa, aa):20s}): {porcentaje:5.2f}% ({conteo_global[aa]:,})")

    # Aminoacido mas y menos frecuente
    mas_frecuente = max(composicion, key=composicion.get)
    menos_frecuente = min(composicion, key=composicion.get)

    print(f"\n  Aminoacido mas frecuente:  {mas_frecuente} ({AA_NOMBRES[mas_frecuente]}) - {composicion[mas_frecuente]}%")
    print(f"  Aminoacido menos frecuente: {menos_frecuente} ({AA_NOMBRES[menos_frecuente]}) - {composicion[menos_frecuente]}%")
    print(f"  Total aminoacidos analizados: {total_aa:,}")

    return {
        "composicion_porcentaje": composicion,
        "conteo_absoluto": dict(conteo_global),
        "total_aminoacidos": total_aa,
        "mas_frecuente": {"aminoacido": mas_frecuente, "nombre": AA_NOMBRES[mas_frecuente], "porcentaje": composicion[mas_frecuente]},
        "menos_frecuente": {"aminoacido": menos_frecuente, "nombre": AA_NOMBRES[menos_frecuente], "porcentaje": composicion[menos_frecuente]}
    }


def calcular_propiedades_fisicoquimicas(proteinas):
    """
    Calcula propiedades fisicoquimicas de cada proteina usando ProteinAnalysis.
    Procesa en batches para gestionar memoria.
    """
    print("\n" + "=" * 70)
    print("PROPIEDADES FISICOQUIMICAS DEL PROTEOMA")
    print("=" * 70)

    propiedades_lista = []
    errores = 0
    batch_size = 200

    for i, p in enumerate(proteinas):
        try:
            seq_str = p["secuencia_proteina"]
            if len(seq_str) < 5:
                continue

            analisis = ProteinAnalysis(seq_str)

            mw = analisis.molecular_weight()
            pi = analisis.isoelectric_point()
            gravy = analisis.gravy()
            aromaticidad = analisis.aromaticity()

            try:
                inestabilidad = analisis.instability_index()
            except Exception:
                inestabilidad = None

            try:
                carga_ph7 = analisis.charge_at_pH(7.0)
            except Exception:
                carga_ph7 = None

            propiedades_lista.append({
                "locus_tag": p["locus_tag"],
                "nombre_gen": p["nombre_gen"],
                "producto": p["producto"],
                "proteina_id": p["proteina_id"],
                "longitud_aa": p["longitud_aa"],
                "peso_molecular_da": round(mw, 2),
                "punto_isoelectrico": round(pi, 2),
                "gravy": round(gravy, 4),
                "aromaticidad": round(aromaticidad, 4),
                "indice_inestabilidad": round(inestabilidad, 2) if inestabilidad is not None else None,
                "carga_ph7": round(carga_ph7, 2) if carga_ph7 is not None else None,
                "es_estable": inestabilidad < 40 if inestabilidad is not None else None,
                "es_hidrofobica": gravy > 0,
                "clasificacion_pi": "acida" if pi < 7 else "basica"
            })

        except Exception as e:
            errores += 1
            if errores <= 3:
                print(f"  [AVISO] Error en {p.get('locus_tag', '?')}: {e}")

        if (i + 1) % batch_size == 0:
            gc.collect()
            print(f"  Procesadas {i + 1:,}/{len(proteinas):,} proteinas...")

    gc.collect()
    print(f"[OK] {len(propiedades_lista):,} proteinas analizadas ({errores} errores)")

    # Estadisticas globales
    longitudes = [p["longitud_aa"] for p in propiedades_lista]
    pesos = [p["peso_molecular_da"] for p in propiedades_lista]
    pis = [p["punto_isoelectrico"] for p in propiedades_lista]
    gravys = [p["gravy"] for p in propiedades_lista]
    inestabilidades = [p["indice_inestabilidad"] for p in propiedades_lista if p["indice_inestabilidad"] is not None]

    proteinas_acidas = sum(1 for p in propiedades_lista if p["clasificacion_pi"] == "acida")
    proteinas_basicas = sum(1 for p in propiedades_lista if p["clasificacion_pi"] == "basica")
    proteinas_hidrofobicas = sum(1 for p in propiedades_lista if p["es_hidrofobica"])
    proteinas_estables = sum(1 for p in propiedades_lista if p.get("es_estable") is True)
    proteinas_inestables = sum(1 for p in propiedades_lista if p.get("es_estable") is False)

    stats = {
        "total_analizadas": len(propiedades_lista),
        "longitud_promedio_aa": round(statistics.mean(longitudes), 1) if longitudes else 0,
        "longitud_mediana_aa": round(statistics.median(longitudes), 1) if longitudes else 0,
        "longitud_min_aa": min(longitudes) if longitudes else 0,
        "longitud_max_aa": max(longitudes) if longitudes else 0,
        "peso_molecular_promedio_da": round(statistics.mean(pesos), 1) if pesos else 0,
        "peso_molecular_mediana_da": round(statistics.median(pesos), 1) if pesos else 0,
        "pi_promedio": round(statistics.mean(pis), 2) if pis else 0,
        "pi_mediana": round(statistics.median(pis), 2) if pis else 0,
        "gravy_promedio": round(statistics.mean(gravys), 4) if gravys else 0,
        "proteinas_acidas": proteinas_acidas,
        "proteinas_basicas": proteinas_basicas,
        "proteinas_hidrofobicas": proteinas_hidrofobicas,
        "proteinas_hidrofilicas": len(propiedades_lista) - proteinas_hidrofobicas,
        "proteinas_estables": proteinas_estables,
        "proteinas_inestables": proteinas_inestables,
        "indice_inestabilidad_promedio": round(statistics.mean(inestabilidades), 2) if inestabilidades else 0
    }

    # Top 10 mas grandes y mas pequenas
    ordenadas = sorted(propiedades_lista, key=lambda x: x["longitud_aa"], reverse=True)
    top_grandes = ordenadas[:10]
    top_pequenas = list(reversed(ordenadas[-10:]))

    print(f"\n  Longitud promedio: {stats['longitud_promedio_aa']} aa")
    print(f"  Peso molecular promedio: {fmt_num(stats['peso_molecular_promedio_da'], 0)} Da")
    print(f"  pI promedio: {stats['pi_promedio']}")
    print(f"  GRAVY promedio: {stats['gravy_promedio']}")
    print(f"  Acidas / Basicas: {proteinas_acidas} / {proteinas_basicas}")
    print(f"  Estables / Inestables: {proteinas_estables} / {proteinas_inestables}")

    return {
        "estadisticas": stats,
        "propiedades_lista": propiedades_lista,
        "top_10_mas_grandes": top_grandes,
        "top_10_mas_pequenas": top_pequenas
    }


def detectar_peptido_senal(proteinas):
    """
    Detecta peptidos senal usando heuristica de region hidrofobica N-terminal.
    """
    print("\n" + "=" * 70)
    print("DETECCION DE PEPTIDOS SENAL")
    print("=" * 70)

    con_senal = []

    for p in proteinas:
        seq = p["secuencia_proteina"]
        if len(seq) < 20:
            continue

        # Criterio: comienza con M y tiene region hidrofobica de >=8 aa en los primeros 30
        if seq[0] != 'M':
            continue

        region_n = seq[:30]
        max_hidrofobico = 0
        actual = 0
        for aa in region_n:
            if aa in AA_HIDROFOBICOS:
                actual += 1
                max_hidrofobico = max(max_hidrofobico, actual)
            else:
                actual = 0

        if max_hidrofobico >= 8:
            con_senal.append({
                "locus_tag": p["locus_tag"],
                "nombre_gen": p["nombre_gen"],
                "producto": p["producto"],
                "longitud_aa": p["longitud_aa"],
                "region_n_terminal": region_n,
                "longitud_hidrofobica": max_hidrofobico
            })

    total = len(proteinas)
    detectadas = len(con_senal)
    porcentaje = (detectadas / total * 100) if total > 0 else 0

    print(f"[OK] {detectadas:,} proteinas con peptido senal ({porcentaje:.1f}%)")
    print(f"  Estas proteinas probablemente son secretadas o de membrana")

    return {
        "total_detectados": detectadas,
        "porcentaje": round(porcentaje, 1),
        "total_analizadas": total,
        "lista": con_senal[:50]  # Top 50 para el JSON
    }


def categorizar_funciones(proteinas):
    """
    Categoriza proteinas segun su funcion anotada en el campo 'producto'.
    """
    print("\n" + "=" * 70)
    print("CATEGORIZACION FUNCIONAL DE PROTEINAS")
    print("=" * 70)

    categorias = {
        "Enzimas": {
            "keywords": ["synthase", "kinase", "reductase", "transferase", "hydrolase",
                         "ligase", "isomerase", "oxidase", "dehydrogenase", "protease",
                         "nuclease", "phosphatase", "catalase", "polymerase", "helicase",
                         "endonuclease", "exonuclease", "methyltransferase", "acetyltransferase"],
            "proteinas": [], "total": 0
        },
        "Transportadores": {
            "keywords": ["transporter", "permease", "porin", "pump", "channel",
                         "symporter", "antiporter", "ABC", "efflux", "import"],
            "proteinas": [], "total": 0
        },
        "Reguladores": {
            "keywords": ["regulator", "repressor", "activator", "transcription factor",
                         "response regulator", "sensor", "anti-sigma", "sigma factor",
                         "two-component"],
            "proteinas": [], "total": 0
        },
        "Estructurales": {
            "keywords": ["ribosomal", "flagellar", "fimbrial", "pilus", "membrane protein",
                         "outer membrane", "inner membrane", "cell wall", "peptidoglycan",
                         "structural"],
            "proteinas": [], "total": 0
        },
        "Virulencia": {
            "keywords": ["toxin", "invasin", "secretion system", "adhesin", "hemolysin",
                         "virulence", "pathogenicity", "type III", "type IV", "type VI",
                         "enterotoxin", "cytolysin"],
            "proteinas": [], "total": 0
        },
        "Chaperones y Folding": {
            "keywords": ["chaperone", "chaperonin", "foldase", "DnaK", "DnaJ", "GroEL",
                         "GroES", "heat shock", "cold shock", "stress"],
            "proteinas": [], "total": 0
        },
        "Replicacion y Reparacion": {
            "keywords": ["DNA polymerase", "DNA repair", "recombinase", "recombination",
                         "mismatch repair", "excision repair", "SOS", "primosome",
                         "topoisomerase", "gyrase"],
            "proteinas": [], "total": 0
        },
        "Hipoteticas": {
            "keywords": ["hypothetical", "uncharacterized", "putative", "predicted protein",
                         "DUF", "unknown function", "conserved protein"],
            "proteinas": [], "total": 0
        }
    }

    sin_categoria = []

    for p in proteinas:
        producto_lower = p["producto"].lower()
        categorizada = False

        for cat_nombre, cat_data in categorias.items():
            for kw in cat_data["keywords"]:
                if kw.lower() in producto_lower:
                    cat_data["total"] += 1
                    if len(cat_data["proteinas"]) < 10:  # Max 10 ejemplos por categoria
                        cat_data["proteinas"].append({
                            "locus_tag": p["locus_tag"],
                            "nombre_gen": p["nombre_gen"],
                            "producto": p["producto"],
                            "longitud_aa": p["longitud_aa"]
                        })
                    categorizada = True
                    break
            if categorizada:
                break

        if not categorizada:
            sin_categoria.append(p)

    # Categoria "Otras"
    categorias["Otras"] = {
        "keywords": [],
        "proteinas": [{"locus_tag": p["locus_tag"], "nombre_gen": p["nombre_gen"],
                       "producto": p["producto"], "longitud_aa": p["longitud_aa"]}
                      for p in sin_categoria[:10]],
        "total": len(sin_categoria)
    }

    total = len(proteinas)
    resultado = {}

    print(f"\n  {'Categoria':<30} {'Total':>8} {'Porcentaje':>10}")
    print("  " + "-" * 50)

    for cat_nombre, cat_data in categorias.items():
        porcentaje = (cat_data["total"] / total * 100) if total > 0 else 0
        resultado[cat_nombre] = {
            "total": cat_data["total"],
            "porcentaje": round(porcentaje, 1),
            "ejemplos": cat_data["proteinas"]
        }
        print(f"  {cat_nombre:<30} {cat_data['total']:>8,} {porcentaje:>9.1f}%")

    return resultado


def analizar_mutaciones_patogenicas(proteinas):
    """
    Busca genes criticos para resistencia a antibioticos y analiza posiciones clave.
    """
    print("\n" + "=" * 70)
    print("ANALISIS DE MUTACIONES PATOGENICAS")
    print("=" * 70)

    # Crear indice por nombre de gen
    indice_genes = {}
    for p in proteinas:
        nombre = p["nombre_gen"].lower()
        if nombre:
            indice_genes[nombre] = p

    resultados = []

    for gen_nombre, gen_info in GENES_CRITICOS.items():
        proteina = indice_genes.get(gen_nombre.lower())

        resultado = {
            "gen": gen_nombre,
            "descripcion": gen_info["descripcion"],
            "encontrado": proteina is not None,
            "longitud_esperada_aa": gen_info["longitud_esperada_aa"],
            "longitud_observada_aa": None,
            "diferencia_longitud": None,
            "posiciones_clave": [],
            "producto": ""
        }

        if proteina:
            resultado["longitud_observada_aa"] = proteina["longitud_aa"]
            resultado["diferencia_longitud"] = proteina["longitud_aa"] - gen_info["longitud_esperada_aa"]
            resultado["producto"] = proteina["producto"]
            resultado["locus_tag"] = proteina["locus_tag"]

            seq = proteina["secuencia_proteina"]
            for pos, descripcion in gen_info["posiciones_clave"].items():
                if pos <= len(seq):
                    aa_actual = seq[pos - 1]  # posiciones 1-based
                    resultado["posiciones_clave"].append({
                        "posicion": pos,
                        "aminoacido_actual": aa_actual,
                        "descripcion_mutacion": descripcion
                    })

            estado = "PRESENTE" if abs(resultado["diferencia_longitud"]) < 20 else "LONGITUD ATIPICA"
            print(f"  {gen_nombre}: {estado} ({proteina['longitud_aa']} aa, esperado {gen_info['longitud_esperada_aa']})")
        else:
            print(f"  {gen_nombre}: NO ENCONTRADO")

        resultados.append(resultado)

    encontrados = sum(1 for r in resultados if r["encontrado"])
    print(f"\n[OK] {encontrados}/{len(GENES_CRITICOS)} genes criticos encontrados")

    return resultados


# =============================================================================
# 2. ESTRUCTURA SECUNDARIA
# =============================================================================

def predecir_estructura_secundaria(proteinas):
    """
    Predice estructura secundaria usando ProteinAnalysis.secondary_structure_fraction().
    """
    print("\n" + "=" * 70)
    print("PREDICCION DE ESTRUCTURA SECUNDARIA")
    print("=" * 70)

    helix_total = 0
    sheet_total = 0
    turn_total = 0
    analisis_lista = []
    errores = 0

    for p in proteinas:
        try:
            seq_str = p["secuencia_proteina"]
            if len(seq_str) < 10:
                continue

            analisis = ProteinAnalysis(seq_str)
            helix, turn, sheet = analisis.secondary_structure_fraction()

            analisis_lista.append({
                "locus_tag": p["locus_tag"],
                "nombre_gen": p["nombre_gen"],
                "producto": p["producto"],
                "longitud_aa": p["longitud_aa"],
                "helix": round(helix * 100, 1),
                "sheet": round(sheet * 100, 1),
                "turn": round(turn * 100, 1)
            })

            helix_total += helix
            sheet_total += sheet
            turn_total += turn

        except Exception:
            errores += 1

    n = len(analisis_lista)
    if n == 0:
        print("[ERROR] No se pudo analizar ninguna proteina")
        return {"promedio_proteoma": {"helix": 0, "sheet": 0, "turn": 0}}

    promedio = {
        "helix": round(helix_total / n * 100, 1),
        "sheet": round(sheet_total / n * 100, 1),
        "turn": round(turn_total / n * 100, 1)
    }

    print(f"[OK] {n:,} proteinas analizadas")
    print(f"  Helice alfa promedio: {promedio['helix']}%")
    print(f"  Lamina beta promedio: {promedio['sheet']}%")
    print(f"  Giros/Coil promedio:  {promedio['turn']}%")

    # Top 10 con mas helix
    top_helix = sorted(analisis_lista, key=lambda x: x["helix"], reverse=True)[:10]
    top_sheet = sorted(analisis_lista, key=lambda x: x["sheet"], reverse=True)[:10]

    print(f"\n  Top proteinas con mas helice-alfa:")
    for t in top_helix[:5]:
        print(f"    {t['nombre_gen'] or t['locus_tag']}: {t['helix']}% - {t['producto']}")

    # Distribucion por rangos
    rangos_helix = {"0-20%": 0, "20-40%": 0, "40-60%": 0, "60-80%": 0, "80-100%": 0}
    for a in analisis_lista:
        h = a["helix"]
        if h < 20: rangos_helix["0-20%"] += 1
        elif h < 40: rangos_helix["20-40%"] += 1
        elif h < 60: rangos_helix["40-60%"] += 1
        elif h < 80: rangos_helix["60-80%"] += 1
        else: rangos_helix["80-100%"] += 1

    return {
        "promedio_proteoma": promedio,
        "total_analizadas": n,
        "distribucion_helix": rangos_helix,
        "top_10_helice": top_helix,
        "top_10_lamina": top_sheet
    }


# =============================================================================
# 3. ESTRUCTURA TERCIARIA
# =============================================================================

def generar_info_terciaria(proteinas):
    """
    Genera informacion sobre estructura terciaria: analisis de cisteinas
    y contenido educativo.
    """
    print("\n" + "=" * 70)
    print("ANALISIS DE ESTRUCTURA TERCIARIA")
    print("=" * 70)

    # Analisis de cisteinas (potenciales puentes disulfuro)
    proteinas_con_cys = []
    for p in proteinas:
        seq = p["secuencia_proteina"]
        num_cys = seq.count('C')
        if num_cys >= 2:
            proteinas_con_cys.append({
                "locus_tag": p["locus_tag"],
                "nombre_gen": p["nombre_gen"],
                "producto": p["producto"],
                "longitud_aa": p["longitud_aa"],
                "num_cisteinas": num_cys,
                "posiciones_cisteina": [i + 1 for i, aa in enumerate(seq) if aa == 'C'],
                "potenciales_puentes": num_cys // 2
            })

    proteinas_con_cys.sort(key=lambda x: x["num_cisteinas"], reverse=True)

    print(f"[OK] {len(proteinas_con_cys):,} proteinas con 2+ cisteinas (potenciales puentes disulfuro)")
    if proteinas_con_cys:
        print(f"  Maximas cisteinas: {proteinas_con_cys[0]['nombre_gen'] or proteinas_con_cys[0]['locus_tag']} ({proteinas_con_cys[0]['num_cisteinas']} Cys)")

    # Contenido educativo
    info_educativa = {
        "titulo": "Estructura Terciaria de Proteinas",
        "descripcion": "La estructura terciaria describe el plegamiento tridimensional completo de una cadena polipeptidica. Es determinada por interacciones entre los grupos R (cadenas laterales) de los aminoacidos.",
        "fuerzas_estabilizadoras": [
            {"nombre": "Interacciones hidrofobicas", "descripcion": "Aminoacidos no polares se agrupan en el interior de la proteina, lejos del agua"},
            {"nombre": "Puentes de hidrogeno", "descripcion": "Entre grupos polares de cadenas laterales"},
            {"nombre": "Puentes disulfuro", "descripcion": "Enlaces covalentes entre cisteinas (-S-S-). Detectamos proteinas con multiples Cys como candidatas"},
            {"nombre": "Interacciones ionicas", "descripcion": "Entre aminoacidos con carga opuesta (Asp/Glu con Lys/Arg)"},
            {"nombre": "Fuerzas de van der Waals", "descripcion": "Interacciones debiles entre atomos cercanos que contribuyen a la estabilidad global"}
        ],
        "metodos_experimentales": [
            {"nombre": "Cristalografia de Rayos X", "descripcion": "Resolucion atomica (~1-3 A). Proteina debe cristalizar. Base de datos: PDB"},
            {"nombre": "Cryo-EM", "descripcion": "Microscopia electronica criogenica. No requiere cristalizacion. Ideal para complejos grandes"},
            {"nombre": "RMN/NMR", "descripcion": "Resonancia magnetica nuclear. Proteina en solucion. Limitado a proteinas pequenas (<40 kDa)"},
            {"nombre": "AlphaFold", "descripcion": "Prediccion computacional con IA (DeepMind). Alta precision. Base de datos: alphafold.ebi.ac.uk"}
        ],
        "nota_bacterias": "En bacterias, la mayoria de proteinas citoplasmaticas no forman puentes disulfuro (ambiente reductor). Las proteinas periplasmaticas y secretadas si pueden tener puentes disulfuro gracias a la maquinaria Dsb (DsbA/DsbB)."
    }

    return {
        "analisis_cisteinas": {
            "total_con_multiples_cys": len(proteinas_con_cys),
            "proteinas": proteinas_con_cys[:30],  # Top 30
            "potenciales_puentes_disulfuro": sum(p["potenciales_puentes"] for p in proteinas_con_cys)
        },
        "info_educativa": info_educativa
    }


def consultar_pdb(proteinas, cache):
    """
    Consulta la base de datos PDB para proteinas relevantes.
    Requiere conexion a internet. Usa cache para evitar consultas repetidas.
    """
    if not HAS_REQUESTS:
        print("  [SKIP] Modulo requests no disponible para consultar PDB")
        return []

    print("\n  Consultando Protein Data Bank (PDB)...")

    # Seleccionar proteinas relevantes (no hipoteticas, las mas largas)
    candidatas = [p for p in proteinas if "hypothetical" not in p["producto"].lower()
                  and "uncharacterized" not in p["producto"].lower()
                  and p["producto"].strip() != ""]
    candidatas.sort(key=lambda x: x["longitud_aa"], reverse=True)
    candidatas = candidatas[:15]  # Solo las 15 mas relevantes

    resultados = []
    pdb_cache = cache.get("pdb", {})

    for p in candidatas:
        locus = p["locus_tag"]

        # Verificar cache
        if locus in pdb_cache:
            resultados.append(pdb_cache[locus])
            continue

        try:
            # Buscar por nombre de proteina
            query = {
                "query": {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": p["producto"]
                    }
                },
                "return_type": "entry",
                "request_options": {
                    "results_content_type": ["experimental"],
                    "return_all_hits": False,
                    "paginate": {"start": 0, "rows": 3}
                }
            }

            resp = requests.post(
                "https://search.rcsb.org/rcsbsearch/v2/query",
                json=query,
                timeout=10
            )

            resultado = {
                "locus_tag": locus,
                "nombre_gen": p["nombre_gen"],
                "producto": p["producto"],
                "pdb_ids": [],
                "encontrado": False
            }

            if resp.status_code == 200:
                data = resp.json()
                hits = data.get("result_set", [])
                if hits:
                    resultado["pdb_ids"] = [h["identifier"] for h in hits[:3]]
                    resultado["encontrado"] = True
                    resultado["total_hits"] = data.get("total_count", 0)

            pdb_cache[locus] = resultado
            resultados.append(resultado)

        except Exception as e:
            pdb_cache[locus] = {
                "locus_tag": locus, "nombre_gen": p["nombre_gen"],
                "producto": p["producto"], "pdb_ids": [], "encontrado": False,
                "error": str(e)
            }

    cache["pdb"] = pdb_cache
    encontradas = sum(1 for r in resultados if r.get("encontrado"))
    print(f"  [OK] {encontradas}/{len(candidatas)} proteinas con estructura PDB conocida")

    return resultados


def consultar_alphafold(proteinas, cache):
    """
    Consulta AlphaFold DB para proteinas con UniProt ID.
    """
    if not HAS_REQUESTS:
        print("  [SKIP] Modulo requests no disponible para consultar AlphaFold")
        return []

    print("\n  Consultando AlphaFold Database...")

    # Buscar proteinas con protein_id que podria mapearse a UniProt
    candidatas = [p for p in proteinas if p["proteina_id"]
                  and "hypothetical" not in p["producto"].lower()][:10]

    resultados = []
    af_cache = cache.get("alphafold", {})

    for p in candidatas:
        locus = p["locus_tag"]

        if locus in af_cache:
            resultados.append(af_cache[locus])
            continue

        try:
            # Intentar buscar por nombre del gen en UniProt
            uniprot_url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{p['nombre_gen']}+AND+organism_id:83333&size=1&format=json"
            resp = requests.get(uniprot_url, timeout=10)

            resultado = {
                "locus_tag": locus,
                "nombre_gen": p["nombre_gen"],
                "producto": p["producto"],
                "uniprot_id": None,
                "alphafold_url": None,
                "encontrado": False
            }

            if resp.status_code == 200:
                data = resp.json()
                results = data.get("results", [])
                if results:
                    uniprot_id = results[0].get("primaryAccession", "")
                    resultado["uniprot_id"] = uniprot_id

                    # Consultar AlphaFold
                    af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
                    af_resp = requests.get(af_url, timeout=10)
                    if af_resp.status_code == 200:
                        af_data = af_resp.json()
                        if af_data:
                            entry = af_data[0] if isinstance(af_data, list) else af_data
                            resultado["alphafold_url"] = entry.get("pdbUrl", "")
                            resultado["confidence"] = entry.get("paeDocUrl", "")
                            resultado["encontrado"] = True

            af_cache[locus] = resultado
            resultados.append(resultado)

        except Exception:
            pass

    cache["alphafold"] = af_cache
    encontradas = sum(1 for r in resultados if r.get("encontrado"))
    print(f"  [OK] {encontradas}/{len(candidatas)} proteinas con modelo AlphaFold")

    return resultados


# =============================================================================
# 4. ESTRUCTURA CUATERNARIA
# =============================================================================

def analizar_estructura_cuaternaria(proteinas):
    """
    Detecta complejos proteicos multi-subunidad analizando anotaciones.
    """
    print("\n" + "=" * 70)
    print("ANALISIS DE ESTRUCTURA CUATERNARIA")
    print("=" * 70)

    # Keywords que indican subunidades
    keywords_subunidad = ["subunit", "subunidad", "component", "chain"]
    keywords_complejo = ["complex", "complejo", "dimer", "trimer", "tetramer",
                         "hexamer", "octamer", "decamer", "dodecamer",
                         "homodimer", "heterodimer", "homotetramer"]
    keywords_tipo = ["alpha", "beta", "gamma", "delta", "epsilon", "omega",
                     "large", "small", "I", "II", "III", "IV"]

    subunidades = []
    complejos_dict = defaultdict(list)

    for p in proteinas:
        producto = p["producto"]
        producto_lower = producto.lower()

        es_subunidad = False
        for kw in keywords_subunidad:
            if kw in producto_lower:
                es_subunidad = True
                break

        if not es_subunidad:
            for kw in keywords_complejo:
                if kw in producto_lower:
                    es_subunidad = True
                    break

        if es_subunidad:
            # Extraer nombre base del complejo
            nombre_base = producto
            # Remover tipo de subunidad
            for tipo in ["subunit alpha", "subunit beta", "subunit gamma", "subunit delta",
                         "subunit epsilon", "subunit omega", "subunit A", "subunit B",
                         "subunit C", "subunit D", "subunit E", "subunit F",
                         "subunit G", "subunit H", "subunit I", "subunit II",
                         "subunit III", "subunit IV", "large subunit", "small subunit",
                         "alpha subunit", "beta subunit", "gamma subunit"]:
                nombre_base = nombre_base.replace(tipo, "").strip()

            # Limpiar nombre
            nombre_base = re.sub(r'\s+', ' ', nombre_base).strip()
            nombre_base = nombre_base.rstrip(',').strip()

            if nombre_base:
                info_sub = {
                    "locus_tag": p["locus_tag"],
                    "nombre_gen": p["nombre_gen"],
                    "producto": producto,
                    "longitud_aa": p["longitud_aa"],
                    "nombre_base_complejo": nombre_base
                }
                subunidades.append(info_sub)
                complejos_dict[nombre_base].append(info_sub)

    # Filtrar complejos con 2+ subunidades
    complejos = []
    for nombre, subs in complejos_dict.items():
        if len(subs) >= 2:
            complejos.append({
                "nombre_complejo": nombre,
                "num_subunidades": len(subs),
                "subunidades": [{
                    "locus_tag": s["locus_tag"],
                    "nombre_gen": s["nombre_gen"],
                    "producto": s["producto"],
                    "longitud_aa": s["longitud_aa"]
                } for s in subs],
                "peso_total_estimado_aa": sum(s["longitud_aa"] for s in subs)
            })

    complejos.sort(key=lambda x: x["num_subunidades"], reverse=True)

    print(f"[OK] {len(complejos)} complejos proteicos detectados")
    print(f"  Total subunidades identificadas: {len(subunidades)}")

    if complejos:
        print(f"\n  Complejos mas grandes:")
        for c in complejos[:10]:
            genes = ", ".join([s["nombre_gen"] or s["locus_tag"] for s in c["subunidades"]])
            print(f"    {c['nombre_complejo']}: {c['num_subunidades']} subunidades ({genes})")

    complejo_mayor = complejos[0] if complejos else None

    return {
        "total_complejos": len(complejos),
        "total_subunidades": len(subunidades),
        "complejos_detectados": complejos[:30],  # Top 30
        "complejo_mas_grande": complejo_mayor
    }


# =============================================================================
# GENERACION DE FIGURAS
# =============================================================================

def generar_figuras(composicion, propiedades, categorias, secundaria, cuaternaria):
    """
    Genera todas las figuras PNG del analisis de proteinas.
    """
    print("\n" + "=" * 70)
    print("GENERANDO FIGURAS")
    print("=" * 70)

    # 1. Distribucion de tamano de proteinas
    try:
        longitudes = [p["longitud_aa"] for p in propiedades["propiedades_lista"]]
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(longitudes, bins=50, color=COLORES['primario'], edgecolor='white', alpha=0.85)
        ax.set_xlabel('Longitud (aminoacidos)', fontsize=11)
        ax.set_ylabel('Numero de proteinas', fontsize=11)
        ax.set_title('Distribucion de Tamano de Proteinas', fontsize=13, fontweight='bold')
        ax.axvline(statistics.mean(longitudes), color='red', linestyle='--', linewidth=1.5,
                   label=f'Promedio: {statistics.mean(longitudes):.0f} aa')
        ax.axvline(statistics.median(longitudes), color='orange', linestyle='--', linewidth=1.5,
                   label=f'Mediana: {statistics.median(longitudes):.0f} aa')
        ax.legend(fontsize=9)
        guardar_figura(fig, 'prot_distribucion_tamano')
    except Exception as e:
        print(f"  [ERROR] Figura tamano: {e}")

    # 2. Composicion de aminoacidos
    try:
        comp = composicion["composicion_porcentaje"]
        aas = sorted(comp.keys())
        valores = [comp[aa] for aa in aas]
        nombres = [f"{aa}\n{AA_NOMBRES.get(aa, '')[:3]}" for aa in aas]

        fig, ax = plt.subplots(figsize=(12, 5))
        barras = ax.bar(nombres, valores, color=[COLORES['primario'] if v > 5 else COLORES['azul'] for v in valores],
                        edgecolor='white')
        ax.set_xlabel('Aminoacido', fontsize=11)
        ax.set_ylabel('Porcentaje (%)', fontsize=11)
        ax.set_title('Composicion de Aminoacidos del Proteoma', fontsize=13, fontweight='bold')

        for bar, val in zip(barras, valores):
            ax.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.15,
                    f'{val:.1f}', ha='center', va='bottom', fontsize=7)

        plt.tight_layout()
        guardar_figura(fig, 'prot_composicion_aminoacidos')
    except Exception as e:
        print(f"  [ERROR] Figura composicion: {e}")

    # 3. Distribucion de peso molecular
    try:
        pesos = [p["peso_molecular_da"] / 1000 for p in propiedades["propiedades_lista"]]  # kDa
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(pesos, bins=50, color=COLORES['secundario'], edgecolor='white', alpha=0.85)
        ax.set_xlabel('Peso Molecular (kDa)', fontsize=11)
        ax.set_ylabel('Numero de proteinas', fontsize=11)
        ax.set_title('Distribucion de Peso Molecular', fontsize=13, fontweight='bold')
        ax.axvline(statistics.mean(pesos), color='red', linestyle='--', linewidth=1.5,
                   label=f'Promedio: {statistics.mean(pesos):.1f} kDa')
        ax.legend(fontsize=9)
        guardar_figura(fig, 'prot_distribucion_peso_molecular')
    except Exception as e:
        print(f"  [ERROR] Figura MW: {e}")

    # 4. Distribucion de punto isoelectrico
    try:
        pis = [p["punto_isoelectrico"] for p in propiedades["propiedades_lista"]]
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.hist(pis, bins=40, color=COLORES['terciario'], edgecolor='white', alpha=0.85)
        ax.set_xlabel('Punto Isoelectrico (pI)', fontsize=11)
        ax.set_ylabel('Numero de proteinas', fontsize=11)
        ax.set_title('Distribucion de Punto Isoelectrico', fontsize=13, fontweight='bold')
        ax.axvline(7.0, color='gray', linestyle=':', linewidth=1, label='pH neutro (7.0)')
        ax.axvline(statistics.mean(pis), color='red', linestyle='--', linewidth=1.5,
                   label=f'Promedio: {statistics.mean(pis):.1f}')
        ax.legend(fontsize=9)
        guardar_figura(fig, 'prot_distribucion_pi')
    except Exception as e:
        print(f"  [ERROR] Figura pI: {e}")

    # 5. Categorias funcionales (donut chart)
    try:
        cats_filtradas = {k: v for k, v in categorias.items() if v["total"] > 0}
        nombres_cat = list(cats_filtradas.keys())
        totales_cat = [cats_filtradas[n]["total"] for n in nombres_cat]

        colores_cat = [COLORES['primario'], COLORES['azul'], COLORES['terciario'],
                       COLORES['secundario'], COLORES['cuaternario'], COLORES['cyan'],
                       COLORES['verde'], COLORES['gris'], '#8B5CF6'][:len(nombres_cat)]

        fig, ax = plt.subplots(figsize=(10, 7))
        wedges, texts, autotexts = ax.pie(
            totales_cat, labels=None, autopct='%1.1f%%',
            colors=colores_cat, startangle=90,
            pctdistance=0.8, wedgeprops=dict(width=0.4, edgecolor='white')
        )

        for t in autotexts:
            t.set_fontsize(7)

        ax.legend(wedges, [f"{n} ({t})" for n, t in zip(nombres_cat, totales_cat)],
                  loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)
        ax.set_title('Categorias Funcionales del Proteoma', fontsize=13, fontweight='bold')
        plt.tight_layout()
        guardar_figura(fig, 'prot_categorias_funcionales')
    except Exception as e:
        print(f"  [ERROR] Figura categorias: {e}")

    # 6. Estructura secundaria del proteoma
    try:
        prom = secundaria["promedio_proteoma"]
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Pie chart de promedios
        valores_sec = [prom["helix"], prom["sheet"], prom["turn"]]
        labels_sec = [f'Helice alfa\n{prom["helix"]}%', f'Lamina beta\n{prom["sheet"]}%', f'Giros/Coil\n{prom["turn"]}%']
        colores_sec = [COLORES['primario'], COLORES['azul'], COLORES['terciario']]

        axes[0].pie(valores_sec, labels=labels_sec, colors=colores_sec,
                    autopct='%1.1f%%', startangle=90, textprops={'fontsize': 9})
        axes[0].set_title('Composicion Promedio\nde Estructura Secundaria', fontsize=11, fontweight='bold')

        # Barras de distribucion de helix
        if "distribucion_helix" in secundaria:
            dist = secundaria["distribucion_helix"]
            rangos = list(dist.keys())
            conteos = list(dist.values())
            axes[1].bar(rangos, conteos, color=COLORES['primario'], edgecolor='white')
            axes[1].set_xlabel('Porcentaje de helice alfa', fontsize=10)
            axes[1].set_ylabel('Numero de proteinas', fontsize=10)
            axes[1].set_title('Distribucion de Contenido\nde Helice Alfa', fontsize=11, fontweight='bold')

        plt.tight_layout()
        guardar_figura(fig, 'prot_estructura_secundaria')
    except Exception as e:
        print(f"  [ERROR] Figura secundaria: {e}")

    # 7. Complejos cuaternarios
    try:
        complejos = cuaternaria.get("complejos_detectados", [])
        if complejos:
            top_complejos = complejos[:15]
            nombres_comp = [c["nombre_complejo"][:30] for c in top_complejos]
            num_subs = [c["num_subunidades"] for c in top_complejos]

            fig, ax = plt.subplots(figsize=(12, 6))
            bars = ax.barh(range(len(nombres_comp)), num_subs, color=COLORES['cuaternario'], edgecolor='white')
            ax.set_yticks(range(len(nombres_comp)))
            ax.set_yticklabels(nombres_comp, fontsize=8)
            ax.set_xlabel('Numero de subunidades', fontsize=11)
            ax.set_title('Complejos Proteicos (por numero de subunidades)', fontsize=13, fontweight='bold')
            ax.invert_yaxis()

            for bar, n in zip(bars, num_subs):
                ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2.,
                        str(n), ha='left', va='center', fontsize=9, fontweight='bold')

            plt.tight_layout()
            guardar_figura(fig, 'prot_complejos_cuaternarios')
        else:
            print("  [SKIP] No hay complejos para graficar")
    except Exception as e:
        print(f"  [ERROR] Figura cuaternaria: {e}")


# =============================================================================
# EXPORTACION DE RESULTADOS
# =============================================================================

def exportar_json(proteinas, composicion, propiedades, peptidos, categorias,
                  mutaciones, secundaria, terciaria, cuaternaria,
                  pdb_resultados, alphafold_resultados, organismo):
    """
    Exporta todos los resultados a un archivo JSON consolidado.
    """
    archivo = os.path.join(RUTA_TABLAS, f"analisis_proteinas_{GENOME_BASENAME}.json")
    print(f"\nExportando JSON: {os.path.basename(archivo)}")

    # Preparar propiedades sin la lista completa (para reducir tamano)
    props_resumen = {
        "estadisticas": propiedades["estadisticas"],
        "top_10_mas_grandes": propiedades["top_10_mas_grandes"],
        "top_10_mas_pequenas": propiedades["top_10_mas_pequenas"]
    }

    datos = {
        "metadata": {
            "fecha_analisis": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "script": "analisis_proteinas.py",
            "genoma": GENOME_BASENAME,
            "organismo": organismo,
            "total_proteinas": len(proteinas)
        },
        "estructura_primaria": {
            "estadisticas_generales": propiedades["estadisticas"],
            "composicion_aminoacidos": composicion,
            "top_10_mas_grandes": propiedades["top_10_mas_grandes"],
            "top_10_mas_pequenas": propiedades["top_10_mas_pequenas"],
            "categorias_funcionales": categorias,
            "peptidos_senal": peptidos,
            "mutaciones_patogenicas": mutaciones
        },
        "estructura_secundaria": secundaria,
        "estructura_terciaria": {
            "analisis_cisteinas": terciaria["analisis_cisteinas"],
            "info_educativa": terciaria["info_educativa"],
            "proteinas_con_pdb": pdb_resultados,
            "proteinas_alphafold": alphafold_resultados
        },
        "estructura_cuaternaria": cuaternaria
    }

    with open(archivo, 'w', encoding='utf-8') as f:
        json.dump(datos, f, ensure_ascii=False, indent=2)

    print(f"[OK] JSON exportado: {archivo}")
    return datos


def exportar_csv_propiedades(propiedades):
    """
    Exporta propiedades fisicoquimicas de todas las proteinas a CSV.
    """
    archivo = os.path.join(RUTA_TABLAS, f"proteinas_propiedades_{GENOME_BASENAME}.csv")
    print(f"Exportando CSV propiedades: {os.path.basename(archivo)}")

    campos = [
        "locus_tag", "nombre_gen", "producto", "proteina_id",
        "longitud_aa", "peso_molecular_da", "punto_isoelectrico",
        "gravy", "aromaticidad", "indice_inestabilidad", "carga_ph7",
        "es_estable", "es_hidrofobica", "clasificacion_pi"
    ]

    with open(archivo, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=campos, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(propiedades["propiedades_lista"])

    print(f"[OK] {len(propiedades['propiedades_lista']):,} proteinas exportadas a CSV")


def exportar_csv_composicion(composicion):
    """
    Exporta composicion de aminoacidos a CSV.
    """
    archivo = os.path.join(RUTA_TABLAS, f"composicion_aminoacidos_{GENOME_BASENAME}.csv")
    print(f"Exportando CSV composicion: {os.path.basename(archivo)}")

    with open(archivo, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["aminoacido", "nombre", "porcentaje", "conteo"])
        for aa in sorted(composicion["composicion_porcentaje"].keys()):
            writer.writerow([
                aa,
                AA_NOMBRES.get(aa, ""),
                composicion["composicion_porcentaje"][aa],
                composicion["conteo_absoluto"].get(aa, 0)
            ])

    print(f"[OK] Composicion de aminoacidos exportada a CSV")


def exportar_csv_complejos(cuaternaria):
    """
    Exporta complejos proteicos a CSV.
    """
    archivo = os.path.join(RUTA_TABLAS, f"complejos_proteicos_{GENOME_BASENAME}.csv")
    print(f"Exportando CSV complejos: {os.path.basename(archivo)}")

    with open(archivo, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["complejo", "num_subunidades", "subunidades_genes", "subunidades_productos", "peso_total_aa"])
        for c in cuaternaria.get("complejos_detectados", []):
            genes = "; ".join([s.get("nombre_gen", "") or s.get("locus_tag", "") for s in c["subunidades"]])
            productos = "; ".join([s.get("producto", "") for s in c["subunidades"]])
            writer.writerow([
                c["nombre_complejo"],
                c["num_subunidades"],
                genes,
                productos,
                c.get("peso_total_estimado_aa", 0)
            ])

    print(f"[OK] {len(cuaternaria.get('complejos_detectados', []))} complejos exportados a CSV")


# =============================================================================
# FUNCION PRINCIPAL
# =============================================================================

def main():
    """
    Ejecuta el pipeline completo de analisis de proteinas.
    """
    print("=" * 70)
    print("  GenomeHub - ANALISIS COMPLETO DE PROTEINAS")
    print(f"  Genoma: {GENOME_BASENAME}")
    print(f"  Archivo: {os.path.basename(archivo_gb)}")
    print(f"  Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 70)

    # Parsear GenBank
    print(f"\nLeyendo archivo GenBank: {os.path.basename(archivo_gb)}")
    registro = SeqIO.read(archivo_gb, "genbank")
    organismo = registro.annotations.get("organism", "Desconocido")
    print(f"  Organismo: {organismo}")
    print(f"  Longitud genoma: {len(registro.seq):,} pb")

    # =============================================
    # 1. EXTRAER Y TRADUCIR PROTEINAS
    # =============================================
    proteinas = extraer_proteinas(registro)

    if not proteinas:
        print("[ERROR] No se encontraron proteinas en el genoma")
        sys.exit(1)

    # =============================================
    # 2. ESTRUCTURA PRIMARIA
    # =============================================
    composicion = analizar_composicion_aminoacidos(proteinas)
    gc.collect()

    propiedades = calcular_propiedades_fisicoquimicas(proteinas)
    gc.collect()

    peptidos = detectar_peptido_senal(proteinas)

    categorias = categorizar_funciones(proteinas)

    mutaciones = analizar_mutaciones_patogenicas(proteinas)

    # =============================================
    # 3. ESTRUCTURA SECUNDARIA
    # =============================================
    secundaria = predecir_estructura_secundaria(proteinas)
    gc.collect()

    # =============================================
    # 4. ESTRUCTURA TERCIARIA
    # =============================================
    terciaria = generar_info_terciaria(proteinas)

    # Consultas externas (opcionales)
    cache = cargar_cache()
    pdb_resultados = []
    alphafold_resultados = []

    try:
        pdb_resultados = consultar_pdb(proteinas, cache)
    except Exception as e:
        print(f"  [AVISO] Error consultando PDB: {e}")

    try:
        alphafold_resultados = consultar_alphafold(proteinas, cache)
    except Exception as e:
        print(f"  [AVISO] Error consultando AlphaFold: {e}")

    guardar_cache(cache)

    # =============================================
    # 5. ESTRUCTURA CUATERNARIA
    # =============================================
    cuaternaria = analizar_estructura_cuaternaria(proteinas)

    # =============================================
    # 6. GENERAR FIGURAS
    # =============================================
    generar_figuras(composicion, propiedades, categorias, secundaria, cuaternaria)

    # =============================================
    # 7. EXPORTAR RESULTADOS
    # =============================================
    exportar_json(proteinas, composicion, propiedades, peptidos, categorias,
                  mutaciones, secundaria, terciaria, cuaternaria,
                  pdb_resultados, alphafold_resultados, organismo)

    exportar_csv_propiedades(propiedades)
    exportar_csv_composicion(composicion)
    exportar_csv_complejos(cuaternaria)

    # =============================================
    # RESUMEN FINAL
    # =============================================
    print("\n" + "=" * 70)
    print("  RESUMEN DEL ANALISIS DE PROTEINAS")
    print("=" * 70)
    stats = propiedades["estadisticas"]
    print(f"  Total proteinas: {stats['total_analizadas']:,}")
    print(f"  Longitud promedio: {stats['longitud_promedio_aa']} aa")
    print(f"  Peso molecular promedio: {fmt_num(stats['peso_molecular_promedio_da'], 0)} Da")
    print(f"  pI promedio: {stats['pi_promedio']}")
    print(f"  Proteinas acidas/basicas: {stats['proteinas_acidas']}/{stats['proteinas_basicas']}")
    print(f"  Estructura secundaria: {secundaria['promedio_proteoma']['helix']}% helix, {secundaria['promedio_proteoma']['sheet']}% sheet")
    print(f"  Complejos cuaternarios: {cuaternaria['total_complejos']}")
    print(f"  Peptidos senal: {peptidos['total_detectados']}")
    print(f"  Genes criticos encontrados: {sum(1 for m in mutaciones if m['encontrado'])}/{len(mutaciones)}")
    print(f"\n  Resultados en: {RUTA_TABLAS}")
    print(f"  Figuras en: {RUTA_FIGURAS}")
    print("=" * 70)
    print("[OK] ANALISIS DE PROTEINAS COMPLETADO EXITOSAMENTE")


if __name__ == "__main__":
    main()

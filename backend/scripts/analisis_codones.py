#!/usr/bin/env python3
"""
Analisis de Codones de Genomas Bacterianos

Este script analiza todos los 64 codones posibles, codones de inicio (ATG)
y parada (TAA, TAG, TGA) en el genoma completo.

Funcionalidades:
- Contar los 64 codones posibles y sus frecuencias
- Contar codones ATG y calcular densidad por kilobase
- Contar codones de parada (TAA, TAG, TGA) y sus proporciones
- Calcular contenido GC
- Exportar resultados en formato CSV y JSON

Uso: python analisis_codones.py <genome_basename>
Ejemplo: python analisis_codones.py ecoli_k12

Fecha: 2026
Proyecto: Analisis comparativo de genomas bacterianos
"""

from Bio import SeqIO
import os
import sys
import json
import csv
from datetime import datetime
from collections import Counter


# CONFIGURACION
# =============================================================================

DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # backend/
DIRECTORIO_PROYECTO = os.path.dirname(DIRECTORIO_BASE)  # raiz del proyecto
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_PROYECTO, "datos", "crudo")
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")

# Parametro de linea de comandos
if len(sys.argv) < 2:
    print("[ERROR] Uso: python analisis_codones.py <genome_basename>")
    print("  Ejemplo: python analisis_codones.py ecoli_k12")
    sys.exit(1)

GENOME_BASENAME = sys.argv[1]
ARCHIVO_FASTA = os.path.join(RUTA_DATOS_CRUDO, f"{GENOME_BASENAME}.fasta")

# Codones
CODON_INICIO = "ATG"
CODONES_PARADA = ["TAA", "TAG", "TGA"]

# Tabla del codigo genetico - los 64 codones y su aminoacido
TABLA_CODONES = {
    "TTT": "Phe (F)", "TTC": "Phe (F)", "TTA": "Leu (L)", "TTG": "Leu (L)",
    "CTT": "Leu (L)", "CTC": "Leu (L)", "CTA": "Leu (L)", "CTG": "Leu (L)",
    "ATT": "Ile (I)", "ATC": "Ile (I)", "ATA": "Ile (I)", "ATG": "Met (M)",
    "GTT": "Val (V)", "GTC": "Val (V)", "GTA": "Val (V)", "GTG": "Val (V)",
    "TCT": "Ser (S)", "TCC": "Ser (S)", "TCA": "Ser (S)", "TCG": "Ser (S)",
    "CCT": "Pro (P)", "CCC": "Pro (P)", "CCA": "Pro (P)", "CCG": "Pro (P)",
    "ACT": "Thr (T)", "ACC": "Thr (T)", "ACA": "Thr (T)", "ACG": "Thr (T)",
    "GCT": "Ala (A)", "GCC": "Ala (A)", "GCA": "Ala (A)", "GCG": "Ala (A)",
    "TAT": "Tyr (Y)", "TAC": "Tyr (Y)", "TAA": "Stop", "TAG": "Stop",
    "CAT": "His (H)", "CAC": "His (H)", "CAA": "Gln (Q)", "CAG": "Gln (Q)",
    "AAT": "Asn (N)", "AAC": "Asn (N)", "AAA": "Lys (K)", "AAG": "Lys (K)",
    "GAT": "Asp (D)", "GAC": "Asp (D)", "GAA": "Glu (E)", "GAG": "Glu (E)",
    "TGT": "Cys (C)", "TGC": "Cys (C)", "TGA": "Stop", "TGG": "Trp (W)",
    "CGT": "Arg (R)", "CGC": "Arg (R)", "CGA": "Arg (R)", "CGG": "Arg (R)",
    "AGT": "Ser (S)", "AGC": "Ser (S)", "AGA": "Arg (R)", "AGG": "Arg (R)",
    "GGT": "Gly (G)", "GGC": "Gly (G)", "GGA": "Gly (G)", "GGG": "Gly (G)",
}

# Valores de referencia (solo para E. coli K-12)
VALORES_LITERATURA_ECOLI = {
    "longitud_genoma": 4641652,
    "genes_anotados": 4300,
    "contenido_gc": 50.8,
    "proporcion_taa": 63.0,
    "proporcion_tag": 8.0,
    "proporcion_tga": 29.0,
}

if "ecoli" in GENOME_BASENAME.lower():
    VALORES_LITERATURA = VALORES_LITERATURA_ECOLI
else:
    VALORES_LITERATURA = {}


# FUNCIONES DE UTILIDAD
# =============================================================================

def crear_directorios():
    """Crea los directorios de resultados si no existen."""
    os.makedirs(RUTA_RESULTADOS, exist_ok=True)


def cargar_secuencia(ruta_archivo):
    """
    Carga la secuencia del genoma desde un archivo FASTA.

    Args:
        ruta_archivo: Ruta al archivo FASTA

    Returns:
        str: Secuencia de ADN en mayusculas
    """
    print(f"[INFO] Cargando secuencia desde: {os.path.basename(ruta_archivo)}")

    if not os.path.exists(ruta_archivo):
        raise FileNotFoundError(f"No se encontro el archivo: {ruta_archivo}")

    registro = SeqIO.read(ruta_archivo, "fasta")
    secuencia = str(registro.seq).upper()

    print(f"[OK] Secuencia cargada: {len(secuencia):,} pares de bases")
    return secuencia


# FUNCIONES DE ANALISIS
# =============================================================================

def contar_codon(secuencia, codon):
    """Cuenta las ocurrencias de un codon en la secuencia (sliding window)."""
    contador = 0
    codon = codon.upper()
    for i in range(len(secuencia) - 2):
        if secuencia[i:i+3] == codon:
            contador += 1
    return contador


def calcular_densidad_por_kilobase(conteo, longitud_genoma):
    """Calcula la densidad de un codon por kilobase de genoma."""
    kilobases = longitud_genoma / 1000
    return conteo / kilobases


def contar_64_codones(secuencia):
    """
    Cuenta los 64 codones en los 3 marcos de lectura de la hebra forward.

    Args:
        secuencia: Secuencia de ADN

    Returns:
        dict: Conteo de cada codon, frecuencias, agrupados por aminoacido
    """
    print("\n" + "=" * 60)
    print("CONTEO DE LOS 64 CODONES")
    print("=" * 60)

    longitud = len(secuencia)

    # Contar en los 3 marcos de lectura (frame 0, 1, 2)
    conteo_total = Counter()

    for frame in range(3):
        for i in range(frame, longitud - 2, 3):
            codon = secuencia[i:i+3]
            if len(codon) == 3 and all(n in "ATGC" for n in codon):
                conteo_total[codon] += 1

    total_codones = sum(conteo_total.values())

    # Construir resultados detallados para los 64 codones
    codones_detalle = {}
    for codon in sorted(TABLA_CODONES.keys()):
        conteo = conteo_total.get(codon, 0)
        frecuencia = (conteo / total_codones * 100) if total_codones > 0 else 0
        codones_detalle[codon] = {
            "conteo": conteo,
            "frecuencia_porcentaje": round(frecuencia, 4),
            "aminoacido": TABLA_CODONES[codon],
            "densidad_por_kb": round(conteo / (longitud / 1000), 4)
        }

    # Agrupar por aminoacido
    por_aminoacido = {}
    for codon, info in codones_detalle.items():
        aa = info["aminoacido"]
        if aa not in por_aminoacido:
            por_aminoacido[aa] = {"codones": {}, "total": 0}
        por_aminoacido[aa]["codones"][codon] = info["conteo"]
        por_aminoacido[aa]["total"] += info["conteo"]

    # Mostrar resumen
    print(f"\n  Total codones contados (3 frames): {total_codones:,}")
    print(f"\n  {'Codon':<8} {'AA':<12} {'Conteo':>10} {'Frecuencia':>12}")
    print("  " + "-" * 45)
    for codon in sorted(TABLA_CODONES.keys()):
        info = codones_detalle[codon]
        print(f"  {codon:<8} {info['aminoacido']:<12} {info['conteo']:>10,} {info['frecuencia_porcentaje']:>11.4f}%")

    resultados = {
        "total_codones_contados": total_codones,
        "codones_detalle": codones_detalle,
        "por_aminoacido": por_aminoacido
    }

    return resultados


def analizar_codones_inicio(secuencia):
    """Analiza los codones de inicio (ATG) en el genoma."""
    print("\n" + "=" * 60)
    print("ANALISIS DE CODONES DE INICIO (ATG)")
    print("=" * 60)

    longitud = len(secuencia)
    conteo_atg = contar_codon(secuencia, CODON_INICIO)
    densidad = calcular_densidad_por_kilobase(conteo_atg, longitud)

    genes_anotados = VALORES_LITERATURA.get("genes_anotados", 0)
    ratio_atg_genes = conteo_atg / genes_anotados if genes_anotados > 0 else 0

    resultados = {
        "codon": CODON_INICIO,
        "conteo_total": conteo_atg,
        "longitud_genoma_pb": longitud,
        "densidad_por_kb": round(densidad, 4),
        "genes_anotados_literatura": genes_anotados,
        "ratio_atg_por_gen": round(ratio_atg_genes, 2),
    }

    print(f"\n  Conteo total de ATG:        {conteo_atg:,}")
    print(f"  Longitud del genoma:        {longitud:,} pb")
    print(f"  Densidad por kilobase:      {densidad:.4f} ATG/kb")
    if genes_anotados > 0:
        print(f"  Genes anotados (literatura): {genes_anotados:,}")
        print(f"  Ratio ATG/gen:              {ratio_atg_genes:.2f}")

    return resultados


def analizar_codones_parada(secuencia):
    """Analiza los codones de parada (TAA, TAG, TGA) en el genoma."""
    print("\n" + "=" * 60)
    print("ANALISIS DE CODONES DE PARADA (TAA, TAG, TGA)")
    print("=" * 60)

    longitud = len(secuencia)
    conteos = {}
    densidades = {}

    for codon in CODONES_PARADA:
        conteos[codon] = contar_codon(secuencia, codon)
        densidades[codon] = calcular_densidad_por_kilobase(conteos[codon], longitud)

    total_parada = sum(conteos.values())
    proporciones = {}
    for codon in CODONES_PARADA:
        proporciones[codon] = (conteos[codon] / total_parada) * 100 if total_parada > 0 else 0

    # Proporciones de literatura si disponibles
    proporciones_lit = {}
    if VALORES_LITERATURA:
        proporciones_lit = {
            "TAA": VALORES_LITERATURA.get("proporcion_taa", 0),
            "TAG": VALORES_LITERATURA.get("proporcion_tag", 0),
            "TGA": VALORES_LITERATURA.get("proporcion_tga", 0),
        }

    resultados = {
        "conteos": conteos,
        "total_codones_parada": total_parada,
        "densidades_por_kb": {k: round(v, 4) for k, v in densidades.items()},
        "proporciones_porcentaje": {k: round(v, 2) for k, v in proporciones.items()},
    }
    if proporciones_lit:
        resultados["proporciones_literatura"] = proporciones_lit

    # Mostrar resultados
    print(f"\n  {'Codon':<8} {'Conteo':>12} {'Densidad/kb':>15} {'Proporcion':>12}")
    print("  " + "-" * 50)
    for codon in CODONES_PARADA:
        print(f"  {codon:<8} {conteos[codon]:>12,} {densidades[codon]:>15.4f} {proporciones[codon]:>11.2f}%")
    print("  " + "-" * 50)
    print(f"  {'TOTAL':<8} {total_parada:>12,}")

    return resultados


def calcular_contenido_gc(secuencia):
    """Calcula el contenido de GC del genoma."""
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

    gc_literatura = VALORES_LITERATURA.get("contenido_gc", 0)
    diferencia = contenido_gc - gc_literatura if gc_literatura > 0 else 0

    resultados = {
        "conteo_nucleotidos": {
            "A": conteo_a,
            "T": conteo_t,
            "G": conteo_g,
            "C": conteo_c
        },
        "contenido_gc_porcentaje": round(contenido_gc, 2),
        "contenido_at_porcentaje": round(contenido_at, 2),
    }

    if gc_literatura > 0:
        resultados["gc_literatura"] = gc_literatura
        resultados["diferencia_con_literatura"] = round(diferencia, 2)

    print(f"\n  Nucleotidos:")
    print(f"    A: {conteo_a:,} ({(conteo_a/longitud)*100:.2f}%)")
    print(f"    T: {conteo_t:,} ({(conteo_t/longitud)*100:.2f}%)")
    print(f"    G: {conteo_g:,} ({(conteo_g/longitud)*100:.2f}%)")
    print(f"    C: {conteo_c:,} ({(conteo_c/longitud)*100:.2f}%)")
    print(f"\n  Contenido GC: {contenido_gc:.2f}%")
    print(f"  Contenido AT: {contenido_at:.2f}%")

    return resultados


# =============================================================================
# FUNCIONES DE ANALISIS AVANZADO
# =============================================================================

# Codones sinonimos agrupados por aminoacido (codigo genetico estandar)
CODONES_POR_AMINOACIDO = {
    "Phe": ["TTT", "TTC"],
    "Leu": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "Ile": ["ATT", "ATC", "ATA"],
    "Met": ["ATG"],
    "Val": ["GTT", "GTC", "GTA", "GTG"],
    "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "Pro": ["CCT", "CCC", "CCA", "CCG"],
    "Thr": ["ACT", "ACC", "ACA", "ACG"],
    "Ala": ["GCT", "GCC", "GCA", "GCG"],
    "Tyr": ["TAT", "TAC"],
    "His": ["CAT", "CAC"],
    "Gln": ["CAA", "CAG"],
    "Asn": ["AAT", "AAC"],
    "Lys": ["AAA", "AAG"],
    "Asp": ["GAT", "GAC"],
    "Glu": ["GAA", "GAG"],
    "Cys": ["TGT", "TGC"],
    "Trp": ["TGG"],
    "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Gly": ["GGT", "GGC", "GGA", "GGG"],
    "Stop": ["TAA", "TAG", "TGA"],
}


def calcular_rscu(conteo_64):
    """
    Calcula el RSCU (Relative Synonymous Codon Usage) para cada codon.

    RSCU = (frecuencia observada del codon) / (frecuencia esperada si no hay sesgo)
    RSCU = 1.0 → sin sesgo (uso equitativo)
    RSCU > 1.0 → codon preferido (sobreusado)
    RSCU < 1.0 → codon evitado (infrausado)

    Returns:
        dict: RSCU por codon + codones preferidos/evitados por aminoacido
    """
    print("\n" + "=" * 60)
    print("RSCU - USO RELATIVO DE CODONES SINONIMOS")
    print("=" * 60)
    print("  (RSCU > 1.0 = codon preferido, RSCU < 1.0 = codon evitado)")

    codones_detalle = conteo_64["codones_detalle"]
    rscu = {}
    preferidos = {}
    evitados = {}

    print(f"\n  {'AA':<5} {'Codon':<8} {'Conteo':>10} {'RSCU':>8} {'Estado':<15}")
    print("  " + "-" * 50)

    for aa, codones in CODONES_POR_AMINOACIDO.items():
        if aa == "Stop" or len(codones) <= 1:
            for c in codones:
                rscu[c] = 1.0
            continue

        # Total para este aminoacido
        total_aa = sum(codones_detalle.get(c, {}).get("conteo", 0) for c in codones)
        n_sinonimos = len(codones)
        esperado = total_aa / n_sinonimos if n_sinonimos > 0 else 0

        mejor_codon = None
        mejor_rscu = 0

        for c in codones:
            observado = codones_detalle.get(c, {}).get("conteo", 0)
            valor_rscu = (observado / esperado) if esperado > 0 else 0
            rscu[c] = round(valor_rscu, 3)

            estado = "PREFERIDO" if valor_rscu > 1.3 else ("evitado" if valor_rscu < 0.7 else "neutral")
            print(f"  {aa:<5} {c:<8} {observado:>10,} {valor_rscu:>8.3f} {estado}")

            if valor_rscu > mejor_rscu:
                mejor_rscu = valor_rscu
                mejor_codon = c

        if mejor_codon:
            preferidos[aa] = {"codon": mejor_codon, "rscu": mejor_rscu}

        # Codon mas evitado
        peor = min(codones, key=lambda c: rscu.get(c, 0))
        if rscu.get(peor, 0) < 0.7:
            evitados[aa] = {"codon": peor, "rscu": rscu.get(peor, 0)}

    # Codones raros (RSCU < 0.5) - problemáticos para expresión heteróloga
    codones_raros = {c: v for c, v in rscu.items() if v < 0.5 and c not in ["TAA", "TAG", "TGA"]}

    print(f"\n  Codones preferidos (RSCU > 1.3): {sum(1 for v in rscu.values() if v > 1.3)}")
    print(f"  Codones evitados (RSCU < 0.7): {sum(1 for v in rscu.values() if v < 0.7)}")
    print(f"  Codones raros (RSCU < 0.5): {len(codones_raros)} - problematicos para expresion heterologa")

    if codones_raros:
        print(f"\n  Codones raros detectados:")
        for c, v in sorted(codones_raros.items(), key=lambda x: x[1]):
            print(f"    {c} ({TABLA_CODONES.get(c, '?')}): RSCU = {v:.3f}")

    return {
        "rscu_por_codon": rscu,
        "codones_preferidos": preferidos,
        "codones_evitados": evitados,
        "codones_raros": codones_raros,
        "total_preferidos": sum(1 for v in rscu.values() if v > 1.3),
        "total_evitados": sum(1 for v in rscu.values() if v < 0.7),
        "total_raros": len(codones_raros),
        "interpretacion": "RSCU mide el sesgo en el uso de codones sinonimos. Valores altos indican codones optimizados para la maquinaria traduccional del organismo."
    }


def calcular_numero_efectivo_codones(conteo_64):
    """
    Calcula el Numero Efectivo de Codones (Nc) del genoma.
    Nc varia de 20 (sesgo maximo, 1 codon por AA) a 61 (sin sesgo, uso equitativo).
    Nc bajo → fuerte sesgo codonico → genes altamente expresados.
    """
    print("\n" + "=" * 60)
    print("NUMERO EFECTIVO DE CODONES (Nc)")
    print("=" * 60)

    codones_detalle = conteo_64["codones_detalle"]

    # Calcular F para cada familia de sinonimos
    # F = sum(pi^2) donde pi = frecuencia relativa del codon i en su familia
    F_valores = []

    for aa, codones in CODONES_POR_AMINOACIDO.items():
        if aa == "Stop" or len(codones) <= 1:
            continue

        total = sum(codones_detalle.get(c, {}).get("conteo", 0) for c in codones)
        if total == 0:
            continue

        pi_squared_sum = 0
        for c in codones:
            n = codones_detalle.get(c, {}).get("conteo", 0)
            pi = n / total
            pi_squared_sum += pi ** 2

        # F homogeneidad para esta familia
        k = len(codones)
        F = (total * pi_squared_sum - 1) / (total - 1) if total > 1 else 1.0 / k
        F_valores.append((aa, k, F))

    # Agrupar por degeneracion
    f_by_deg = {}
    for aa, k, f in F_valores:
        if k not in f_by_deg:
            f_by_deg[k] = []
        f_by_deg[k].append(f)

    # Calcular Nc segun Wright (1990)
    # Nc = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6 (para familias de 2,3,4,6 codones)
    nc = 2.0  # Met + Trp (1 codon cada uno)

    for deg in [2, 3, 4, 6]:
        if deg in f_by_deg and f_by_deg[deg]:
            f_promedio = statistics.mean(f_by_deg[deg])
            n_familias = {2: 9, 3: 1, 4: 5, 6: 3}.get(deg, 0)
            nc += n_familias / f_promedio if f_promedio > 0 else n_familias

    nc = min(nc, 61)  # Maximo teorico

    interpretacion = ""
    if nc < 35:
        interpretacion = "Sesgo codonico FUERTE - organismo con codones muy optimizados para alta expresion"
    elif nc < 45:
        interpretacion = "Sesgo codonico MODERADO - optimizacion parcial de codones"
    else:
        interpretacion = "Sesgo codonico DEBIL - uso relativamente uniforme de codones sinonimos"

    print(f"\n  Nc = {nc:.1f} (rango: 20 = maximo sesgo, 61 = sin sesgo)")
    print(f"  {interpretacion}")

    return {
        "nc": round(nc, 1),
        "rango": "20 (maximo sesgo) - 61 (sin sesgo)",
        "interpretacion": interpretacion,
        "referencia": "Wright F (1990) Gene 87:23-29"
    }


def analizar_sesgo_por_aminoacido(conteo_64):
    """
    Para cada aminoacido con codones sinonimos, muestra cual es el preferido
    y la relacion con el contenido GC del organismo.
    """
    print("\n" + "=" * 60)
    print("SESGO CODONICO POR AMINOACIDO")
    print("=" * 60)

    codones_detalle = conteo_64["codones_detalle"]
    sesgo = {}

    for aa, codones in CODONES_POR_AMINOACIDO.items():
        if aa == "Stop" or len(codones) <= 1:
            continue

        total = sum(codones_detalle.get(c, {}).get("conteo", 0) for c in codones)
        if total == 0:
            continue

        detalles = []
        for c in codones:
            n = codones_detalle.get(c, {}).get("conteo", 0)
            pct = (n / total * 100) if total > 0 else 0
            # GC en 3ra posicion
            gc3 = "GC" if c[2] in "GC" else "AT"
            detalles.append({
                "codon": c,
                "conteo": n,
                "porcentaje": round(pct, 1),
                "termina_gc": gc3 == "GC"
            })

        detalles.sort(key=lambda x: x["conteo"], reverse=True)
        preferido = detalles[0]
        evitado = detalles[-1]

        sesgo[aa] = {
            "codones": detalles,
            "preferido": preferido["codon"],
            "preferido_porcentaje": preferido["porcentaje"],
            "evitado": evitado["codon"],
            "evitado_porcentaje": evitado["porcentaje"],
            "sesgo_gc3": sum(1 for d in detalles if d["termina_gc"] and d["porcentaje"] > 100/len(codones))
        }

    print(f"\n  {'AA':<5} {'Preferido':<12} {'%':>6} {'Evitado':<12} {'%':>6}")
    print("  " + "-" * 45)
    for aa in sorted(sesgo.keys()):
        s = sesgo[aa]
        print(f"  {aa:<5} {s['preferido']:<12} {s['preferido_porcentaje']:>5.1f}% {s['evitado']:<12} {s['evitado_porcentaje']:>5.1f}%")

    return sesgo


# =============================================================================
# FUNCIONES DE EXPORTACION
# =============================================================================

def exportar_json(datos, nombre_archivo):
    """Exporta los resultados a formato JSON."""
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.json")
    with open(ruta, 'w', encoding='utf-8') as archivo:
        json.dump(datos, archivo, indent=2, ensure_ascii=False)
    print(f"[OK] Exportado: {nombre_archivo}.json")


def exportar_csv_codones_64(conteo_64, nombre_archivo):
    """Exporta el conteo de los 64 codones a CSV."""
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")
    with open(ruta, 'w', newline='', encoding='utf-8') as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow(['Codon', 'Aminoacido', 'Conteo', 'Frecuencia_porcentaje', 'Densidad_por_kb'])
        for codon in sorted(TABLA_CODONES.keys()):
            info = conteo_64["codones_detalle"][codon]
            escritor.writerow([
                codon,
                info["aminoacido"],
                info["conteo"],
                info["frecuencia_porcentaje"],
                info["densidad_por_kb"]
            ])
    print(f"[OK] Exportado: {nombre_archivo}.csv")


def exportar_csv_codones(resultados_inicio, resultados_parada, nombre_archivo):
    """Exporta resultados de codones inicio/parada a CSV."""
    ruta = os.path.join(RUTA_RESULTADOS, f"{nombre_archivo}.csv")
    with open(ruta, 'w', newline='', encoding='utf-8') as archivo:
        escritor = csv.writer(archivo)
        escritor.writerow(['Tipo', 'Codon', 'Conteo', 'Densidad_por_kb', 'Proporcion_porcentaje'])
        escritor.writerow([
            'Inicio', resultados_inicio['codon'],
            resultados_inicio['conteo_total'],
            resultados_inicio['densidad_por_kb'], '100.00'
        ])
        for codon in CODONES_PARADA:
            escritor.writerow([
                'Parada', codon,
                resultados_parada['conteos'][codon],
                resultados_parada['densidades_por_kb'][codon],
                resultados_parada['proporciones_porcentaje'][codon]
            ])
    print(f"[OK] Exportado: {nombre_archivo}.csv")


# FUNCION PRINCIPAL
# =============================================================================

def main():
    """Funcion principal que ejecuta todo el analisis de codones."""

    print("\n" + "=" * 60)
    print(f"ANALISIS DE CODONES - {GENOME_BASENAME}")
    print("=" * 60)
    print(f"Fecha de analisis: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    crear_directorios()

    try:
        secuencia = cargar_secuencia(ARCHIVO_FASTA)
    except FileNotFoundError as error:
        print(f"[ERROR] {error}")
        print("[INFO] Ejecuta primero descargar_genoma.py para obtener el genoma")
        return

    # Analisis de los 64 codones
    conteo_64 = contar_64_codones(secuencia)

    # Analisis especificos
    resultados_inicio = analizar_codones_inicio(secuencia)
    resultados_parada = analizar_codones_parada(secuencia)
    resultados_gc = calcular_contenido_gc(secuencia)

    # Analisis avanzado: RSCU
    rscu = calcular_rscu(conteo_64)

    # Numero efectivo de codones
    nc = calcular_numero_efectivo_codones(conteo_64)

    # Sesgo por aminoacido
    sesgo = analizar_sesgo_por_aminoacido(conteo_64)

    # Compilar todos los resultados
    todos_resultados = {
        "fecha_analisis": datetime.now().isoformat(),
        "genoma": GENOME_BASENAME,
        "archivo_fuente": os.path.basename(ARCHIVO_FASTA),
        "conteo_64_codones": conteo_64,
        "codones_inicio": resultados_inicio,
        "codones_parada": resultados_parada,
        "contenido_gc": resultados_gc,
        "rscu": rscu,
        "numero_efectivo_codones": nc,
        "sesgo_por_aminoacido": sesgo,
    }
    if VALORES_LITERATURA:
        todos_resultados["valores_literatura"] = VALORES_LITERATURA

    # Exportar resultados con sufijo del genoma
    print("\n" + "=" * 60)
    print("EXPORTANDO RESULTADOS")
    print("=" * 60)

    exportar_json(todos_resultados, f"analisis_codones_{GENOME_BASENAME}")
    exportar_csv_codones_64(conteo_64, f"codones_64_{GENOME_BASENAME}")
    exportar_csv_codones(resultados_inicio, resultados_parada, f"codones_conteo_{GENOME_BASENAME}")

    # Resumen final
    print("\n" + "=" * 60)
    print("RESUMEN FINAL")
    print("=" * 60)
    print(f"  Genoma analizado:     {GENOME_BASENAME}")
    print(f"  Longitud:             {len(secuencia):,} pb")
    print(f"  Contenido GC:         {resultados_gc['contenido_gc_porcentaje']}%")
    print(f"  Total ATG:            {resultados_inicio['conteo_total']:,}")
    print(f"  Total codones parada: {resultados_parada['total_codones_parada']:,}")
    print(f"  Codones unicos:       {len([c for c in conteo_64['codones_detalle'] if conteo_64['codones_detalle'][c]['conteo'] > 0])}/64")
    print("=" * 60)

    print("\n[OK] Analisis completado exitosamente\n")


if __name__ == "__main__":
    main()

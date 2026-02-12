#!/usr/bin/env python3
"""
ANALISIS DE EVOLUCION GENOMICA / PANGENOMA
===========================================
Compara multiples genomas para identificar:
- Pan-genoma (todos los genes unicos)
- Core genome (genes compartidos por todos)
- Accessory genome (genes en algunos pero no todos)
- Genes unicos por genoma
- Distancias Jaccard entre genomas
- Arbol UPGMA (clustering jerarquico)
- Curva del pan-genoma

Uso: python analisis_evolucion.py genome1,genome2,genome3
     python analisis_evolucion.py all
"""

from Bio import SeqIO
import os
import sys
import json
import re
import math
from datetime import datetime
from collections import defaultdict

# =============================================================================
# CONFIGURACION
# =============================================================================

DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # backend/
DIRECTORIO_PROYECTO = os.path.dirname(DIRECTORIO_BASE)  # raiz del proyecto
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_PROYECTO, "datos", "crudo")
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "evolucion")

os.makedirs(RUTA_RESULTADOS, exist_ok=True)

if len(sys.argv) < 2:
    print("[ERROR] Uso: python analisis_evolucion.py <genomas separados por coma o 'all'>")
    sys.exit(1)

arg_genomas = sys.argv[1]


def obtener_genomas_disponibles():
    """Lista todos los archivos .gb disponibles."""
    genomas = []
    if os.path.exists(RUTA_DATOS_CRUDO):
        for f in os.listdir(RUTA_DATOS_CRUDO):
            if f.endswith(".gb"):
                genomas.append(f.replace(".gb", ""))
    return sorted(genomas)


def normalizar_producto(producto):
    """Normaliza nombre de producto para matching."""
    if not producto:
        return ""
    p = producto.lower().strip()
    # Remover prefijos comunes que no aportan a la identidad
    for prefix in ["putative ", "probable ", "predicted ", "uncharacterized "]:
        if p.startswith(prefix):
            p = p[len(prefix):]
    # Remover sufijos comunes
    p = re.sub(r'\s*\(fragment\)\s*$', '', p)
    p = re.sub(r'\s*-like\s*$', '', p)
    return p.strip()


def es_hypothetical(producto):
    """Detecta si un producto es hypothetical/desconocido."""
    if not producto:
        return True
    p = producto.lower()
    return any(term in p for term in [
        "hypothetical protein", "unknown function", "uncharacterized protein",
        "domain-containing protein", "duf", "predicted protein"
    ])


def extraer_genes_genoma(basename):
    """Extrae CDS de un genoma GenBank."""
    archivo = os.path.join(RUTA_DATOS_CRUDO, f"{basename}.gb")
    if not os.path.exists(archivo):
        print(f"  [WARN] Archivo no encontrado: {archivo}")
        return None

    try:
        registro = SeqIO.read(archivo, "genbank")
    except Exception as e:
        print(f"  [WARN] Error leyendo {basename}: {e}")
        return None

    nombre = registro.description or basename
    longitud = len(registro.seq)
    gc_count = sum(1 for n in str(registro.seq).upper() if n in "GC")
    gc_pct = round(gc_count / longitud * 100, 2) if longitud > 0 else 0

    genes = []
    for feature in registro.features:
        if feature.type == "CDS":
            producto = feature.qualifiers.get("product", [""])[0]
            nombre_gen = feature.qualifiers.get("gene", [""])[0]
            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
            inicio = int(feature.location.start)
            fin = int(feature.location.end)
            longitud_pb = fin - inicio

            # Contenido GC del gen
            try:
                seq_gen = str(feature.location.extract(registro.seq)).upper()
                gc_gen = sum(1 for n in seq_gen if n in "GC")
                gc_gen_pct = round(gc_gen / len(seq_gen) * 100, 2) if seq_gen else 0
            except Exception:
                gc_gen_pct = 0

            genes.append({
                "producto": producto,
                "producto_normalizado": normalizar_producto(producto),
                "nombre_gen": nombre_gen,
                "locus_tag": locus_tag,
                "inicio": inicio,
                "fin": fin,
                "longitud_pb": longitud_pb,
                "contenido_gc": gc_gen_pct,
                "es_hypothetical": es_hypothetical(producto)
            })

    return {
        "basename": basename,
        "nombre": nombre,
        "total_genes": len(genes),
        "longitud_pb": longitud,
        "gc_porcentaje": gc_pct,
        "genes": genes
    }


def construir_pangenoma(datos_genomas):
    """Construye el pan-genoma a partir de multiples genomas."""
    # Diccionario global: producto_normalizado → set de basenames
    producto_a_genomas = defaultdict(set)
    # Para hypothetical: usar locus_tag como ID unico
    hyp_a_genomas = defaultdict(set)

    # Detalles por genoma
    genes_por_genoma = {}

    for gdata in datos_genomas:
        basename = gdata["basename"]
        genes_por_genoma[basename] = {
            "total": gdata["total_genes"],
            "productos_set": set()
        }

        for gen in gdata["genes"]:
            if gen["es_hypothetical"]:
                # Los hypothetical se cuentan pero no se usan para matching
                key = f"hyp_{basename}_{gen['locus_tag']}"
                hyp_a_genomas[key].add(basename)
            else:
                prod_norm = gen["producto_normalizado"]
                if prod_norm:
                    producto_a_genomas[prod_norm].add(basename)
                    genes_por_genoma[basename]["productos_set"].add(prod_norm)

    n_genomas = len(datos_genomas)
    basenames = [g["basename"] for g in datos_genomas]

    # Clasificar genes
    core_genes = []
    accessory_genes = []
    unique_genes_map = defaultdict(list)

    for prod, genomas_set in producto_a_genomas.items():
        if len(genomas_set) == n_genomas:
            core_genes.append(prod)
        elif len(genomas_set) == 1:
            genome_owner = list(genomas_set)[0]
            unique_genes_map[genome_owner].append(prod)
        else:
            accessory_genes.append(prod)

    # Hypothetical proteins son siempre unicas (no se matchean)
    total_hyp = len(hyp_a_genomas)

    # Genes unicos detalle
    genes_unicos_detalle = {}
    for gdata in datos_genomas:
        basename = gdata["basename"]
        unicos = unique_genes_map.get(basename, [])
        detalle = []
        for gen in gdata["genes"]:
            if not gen["es_hypothetical"] and gen["producto_normalizado"] in unicos:
                detalle.append({
                    "producto": gen["producto"],
                    "nombre_gen": gen["nombre_gen"],
                    "locus_tag": gen["locus_tag"],
                    "longitud_pb": gen["longitud_pb"]
                })
        genes_unicos_detalle[basename] = detalle

    # Stats por genoma
    stats_por_genoma = {}
    for gdata in datos_genomas:
        basename = gdata["basename"]
        n_unicos = len(unique_genes_map.get(basename, []))
        n_hyp = sum(1 for g in gdata["genes"] if g["es_hypothetical"])
        n_core = sum(1 for g in gdata["genes"]
                     if not g["es_hypothetical"]
                     and g["producto_normalizado"] in [c for c in core_genes])
        stats_por_genoma[basename] = {
            "total": gdata["total_genes"],
            "unicos": n_unicos,
            "core": n_core,
            "hypothetical": n_hyp,
            "porcentaje_core": round(n_core / gdata["total_genes"] * 100, 1) if gdata["total_genes"] > 0 else 0
        }

    total_productos = len(producto_a_genomas)

    return {
        "total_genes_unicos": total_productos + total_hyp,
        "total_productos_conocidos": total_productos,
        "core_genome": len(core_genes),
        "accessory_genome": len(accessory_genes),
        "genes_unicos_total": sum(len(v) for v in unique_genes_map.values()),
        "hypothetical_total": total_hyp,
        "core_genes": sorted(core_genes)[:200],  # Limitar para JSON
        "genes_por_genoma": stats_por_genoma,
        "genes_unicos_detalle": {k: v[:50] for k, v in genes_unicos_detalle.items()},
        "producto_a_genomas": producto_a_genomas,  # temporal, no se exporta
    }


def construir_matriz_presencia(datos_genomas, producto_a_genomas):
    """Construye la matriz de presencia/ausencia."""
    basenames = [g["basename"] for g in datos_genomas]
    # Solo top 500 productos mas frecuentes para no sobrecargar JSON
    productos = sorted(producto_a_genomas.keys(),
                       key=lambda p: len(producto_a_genomas[p]),
                       reverse=True)[:500]

    matriz = []
    for prod in productos:
        fila = [1 if bn in producto_a_genomas[prod] else 0 for bn in basenames]
        matriz.append(fila)

    return {
        "genomas": basenames,
        "genes": productos,
        "matriz": matriz
    }


def calcular_distancias_jaccard(datos_genomas, producto_a_genomas):
    """Calcula matriz de distancias Jaccard entre genomas."""
    basenames = [g["basename"] for g in datos_genomas]

    # Conjuntos de productos por genoma
    sets = {}
    for gdata in datos_genomas:
        bn = gdata["basename"]
        sets[bn] = set()
        for gen in gdata["genes"]:
            if not gen["es_hypothetical"] and gen["producto_normalizado"]:
                sets[bn].add(gen["producto_normalizado"])

    n = len(basenames)
    matriz = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            si = sets[basenames[i]]
            sj = sets[basenames[j]]
            interseccion = len(si & sj)
            union = len(si | sj)
            distancia = 1.0 - (interseccion / union) if union > 0 else 1.0
            matriz[i][j] = round(distancia, 4)
            matriz[j][i] = round(distancia, 4)

    return {
        "genomas": basenames,
        "matriz": matriz
    }


def construir_arbol_upgma(distancias):
    """Construye arbol UPGMA a partir de matriz de distancias."""
    genomas = distancias["genomas"][:]
    mat = [row[:] for row in distancias["matriz"]]
    n = len(genomas)

    # Inicializar nodos hoja
    nodos = []
    for i, nombre in enumerate(genomas):
        nodos.append({
            "id": i,
            "nombre": nombre,
            "hoja": True,
            "hijos": [],
            "distancia": 0,
            "tamano": 1
        })

    # Clusters activos
    activos = list(range(n))
    next_id = n

    while len(activos) > 1:
        # Encontrar par con minima distancia
        min_dist = float('inf')
        mi, mj = 0, 1
        for i in range(len(activos)):
            for j in range(i + 1, len(activos)):
                d = mat[activos[i]][activos[j]]
                if d < min_dist:
                    min_dist = d
                    mi, mj = i, j

        ci = activos[mi]
        cj = activos[mj]

        # Crear nodo interno
        nuevo_nodo = {
            "id": next_id,
            "nombre": None,
            "hoja": False,
            "hijos": [ci, cj],
            "distancia": round(min_dist / 2, 4),
            "tamano": nodos[ci]["tamano"] + nodos[cj]["tamano"]
        }
        nodos.append(nuevo_nodo)

        # Actualizar distancias (UPGMA: promedio ponderado)
        new_row = [0.0] * (next_id + 1)
        for k in activos:
            if k != ci and k != cj:
                ni = nodos[ci]["tamano"]
                nj = nodos[cj]["tamano"]
                d = (mat[ci][k] * ni + mat[cj][k] * nj) / (ni + nj)
                new_row[k] = d
                # Expandir mat si necesario
                while len(mat[k]) <= next_id:
                    mat[k].append(0.0)
                mat[k][next_id] = d

        # Asegurar que mat tiene suficientes filas
        while len(mat) <= next_id:
            mat.append([0.0] * (next_id + 1))
        mat[next_id] = new_row

        # Actualizar activos
        activos.remove(ci)
        activos.remove(cj)
        activos.append(next_id)
        next_id += 1

    # Generar Newick string
    def to_newick(nid):
        n = nodos[nid]
        if n["hoja"]:
            return n["nombre"]
        left = to_newick(n["hijos"][0])
        right = to_newick(n["hijos"][1])
        dl = round(n["distancia"] - nodos[n["hijos"][0]].get("distancia", 0), 4)
        dr = round(n["distancia"] - nodos[n["hijos"][1]].get("distancia", 0), 4)
        return f"({left}:{dl},{right}:{dr})"

    raiz = activos[0] if activos else 0
    newick = to_newick(raiz) + ";"

    # Serializar nodos (sin tamano interno)
    nodos_export = []
    for n in nodos:
        nodos_export.append({
            "id": n["id"],
            "nombre": n["nombre"],
            "hoja": n["hoja"],
            "hijos": n["hijos"],
            "distancia": n["distancia"]
        })

    return {
        "newick": newick,
        "nodos": nodos_export,
        "raiz": raiz
    }


def calcular_curva_pangenoma(datos_genomas):
    """Calcula curva del pan-genoma y core genome."""
    import random
    random.seed(42)

    curva = []
    basenames = [g["basename"] for g in datos_genomas]
    genes_map = {}  # basename → set de productos
    for gdata in datos_genomas:
        s = set()
        for gen in gdata["genes"]:
            if not gen["es_hypothetical"] and gen["producto_normalizado"]:
                s.add(gen["producto_normalizado"])
        genes_map[gdata["basename"]] = s

    # Promediar sobre varias permutaciones
    n_perms = min(20, math.factorial(min(len(basenames), 6)))
    pan_acum = defaultdict(list)
    core_acum = defaultdict(list)

    for _ in range(n_perms):
        orden = basenames[:]
        random.shuffle(orden)
        pan_set = set()
        core_set = None

        for i, bn in enumerate(orden):
            pan_set = pan_set | genes_map[bn]
            if core_set is None:
                core_set = genes_map[bn].copy()
            else:
                core_set = core_set & genes_map[bn]
            pan_acum[i + 1].append(len(pan_set))
            core_acum[i + 1].append(len(core_set))

    for n in range(1, len(basenames) + 1):
        curva.append({
            "n_genomas": n,
            "pan": round(sum(pan_acum[n]) / len(pan_acum[n])),
            "core": round(sum(core_acum[n]) / len(core_acum[n]))
        })

    return curva


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("ANALISIS DE EVOLUCION GENOMICA / PANGENOMA")
    print("=" * 70)
    print(f"  Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Determinar genomas a analizar
    if arg_genomas.lower() == "all":
        genomas_lista = obtener_genomas_disponibles()
    else:
        genomas_lista = [g.strip() for g in arg_genomas.split(",") if g.strip()]

    if len(genomas_lista) < 2:
        print("[ERROR] Se necesitan al menos 2 genomas para el analisis de evolucion")
        sys.exit(1)

    print(f"  Genomas a analizar: {len(genomas_lista)}")
    for g in genomas_lista:
        print(f"    - {g}")

    # Paso 1: Extraer genes de cada genoma
    print(f"\n[PASO 1/{5}] Extrayendo genes de cada genoma...")
    datos_genomas = []
    for i, basename in enumerate(genomas_lista):
        print(f"  [{i+1}/{len(genomas_lista)}] {basename}...")
        gdata = extraer_genes_genoma(basename)
        if gdata:
            datos_genomas.append(gdata)
            print(f"    OK: {gdata['total_genes']} genes, {gdata['longitud_pb']:,} pb, GC={gdata['gc_porcentaje']}%")
        else:
            print(f"    OMITIDO: no se pudo leer")

    if len(datos_genomas) < 2:
        print("[ERROR] Se necesitan al menos 2 genomas validos")
        sys.exit(1)

    # Paso 2: Construir pan-genoma
    print(f"\n[PASO 2/{5}] Construyendo pan-genoma...")
    pangenoma = construir_pangenoma(datos_genomas)
    producto_a_genomas = pangenoma.pop("producto_a_genomas")
    print(f"  Pan-genoma total: {pangenoma['total_genes_unicos']} genes unicos")
    print(f"  Core genome: {pangenoma['core_genome']} genes (presentes en TODOS)")
    print(f"  Accessory: {pangenoma['accessory_genome']} genes (en algunos)")
    print(f"  Unicos: {pangenoma['genes_unicos_total']} genes (en solo 1 genoma)")
    print(f"  Hypothetical: {pangenoma['hypothetical_total']} (no matcheados)")

    # Paso 3: Matriz de presencia/ausencia
    print(f"\n[PASO 3/{5}] Construyendo matriz de presencia/ausencia...")
    matriz_presencia = construir_matriz_presencia(datos_genomas, producto_a_genomas)
    print(f"  Matriz: {len(matriz_presencia['genes'])} genes x {len(matriz_presencia['genomas'])} genomas")

    # Paso 4: Distancias Jaccard + Arbol UPGMA
    print(f"\n[PASO 4/{5}] Calculando distancias Jaccard y arbol UPGMA...")
    distancias = calcular_distancias_jaccard(datos_genomas, producto_a_genomas)

    # Mostrar distancias
    for i in range(len(distancias["genomas"])):
        for j in range(i + 1, len(distancias["genomas"])):
            print(f"  {distancias['genomas'][i][:30]} <-> {distancias['genomas'][j][:30]}: {distancias['matriz'][i][j]:.4f}")

    arbol = construir_arbol_upgma(distancias)
    print(f"  Arbol UPGMA: {len(arbol['nodos'])} nodos")
    print(f"  Newick: {arbol['newick'][:100]}...")

    # Paso 5: Curva del pan-genoma
    print(f"\n[PASO 5/{5}] Calculando curva del pan-genoma...")
    curva = calcular_curva_pangenoma(datos_genomas)
    for punto in curva:
        print(f"  {punto['n_genomas']} genomas: pan={punto['pan']}, core={punto['core']}")

    # Exportar resultado
    resultado = {
        "fecha_analisis": datetime.now().isoformat(),
        "total_genomas": len(datos_genomas),
        "genomas_analizados": [
            {
                "basename": g["basename"],
                "nombre": g["nombre"][:80],
                "total_genes": g["total_genes"],
                "longitud_pb": g["longitud_pb"],
                "gc_porcentaje": g["gc_porcentaje"]
            }
            for g in datos_genomas
        ],
        "pangenoma": pangenoma,
        "matriz_presencia": matriz_presencia,
        "distancias_jaccard": distancias,
        "arbol_upgma": arbol,
        "curva_pangenoma": curva
    }

    archivo_salida = os.path.join(RUTA_RESULTADOS, "analisis_evolucion.json")
    with open(archivo_salida, "w", encoding="utf-8") as f:
        json.dump(resultado, f, ensure_ascii=False, indent=2)

    print(f"\n[OK] Resultado guardado en: {archivo_salida}")
    print(f"  Tamano: {os.path.getsize(archivo_salida) / 1024:.1f} KB")
    print("=" * 70)
    print("[OK] ANALISIS DE EVOLUCION COMPLETADO")
    print("=" * 70)


if __name__ == "__main__":
    main()

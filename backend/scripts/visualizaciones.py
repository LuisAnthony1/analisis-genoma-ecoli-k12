#!/usr/bin/env python3
"""
Visualizaciones del Analisis de Genomas Bacterianos

Este script genera graficos profesionales para ilustrar los resultados
del analisis de genomas bacterianos, incluyendo comparaciones entre
E. coli K-12 MG1655 y Salmonella enterica LT2.

Graficos generados:
- Distribucion de tamanos de genes (histograma)
- Proporciones de codones de parada (grafico circular y barras)
- Contenido GC por gen (histograma)
- Distribucion de genes por hebra (barras)
- Comparacion de resultados con literatura cientifica
- NUEVO: Comparacion E. coli vs Salmonella (metricas, virulencia, distancias)

Fecha: 2026
Proyecto: Analisis comparativo de genomas bacterianos
"""

import matplotlib
matplotlib.use('Agg')  # Backend sin GUI (necesario para AWS)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import json
import os
import sys
import gc
from datetime import datetime


# CONFIGURACION
# =============================================================================

# Rutas de archivos
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # backend/
RUTA_RESULTADOS_TABLAS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")
RUTA_RESULTADOS_FIGURAS = os.path.join(DIRECTORIO_BASE, "resultados", "figuras")

# Parametro de linea de comandos (basename del genoma)
GENOME_BASENAME = sys.argv[1] if len(sys.argv) >= 2 else None

# Archivos de entrada - buscar en carpeta del genoma primero, luego en tablas legacy
def buscar_archivo(nombre_archivo):
    """Busca archivo en carpeta del genoma o en tablas legacy."""
    if GENOME_BASENAME:
        # Primero buscar en resultados/{genome}/tablas/
        ruta_genoma = os.path.join(DIRECTORIO_BASE, "resultados", GENOME_BASENAME, "tablas", nombre_archivo)
        if os.path.exists(ruta_genoma):
            return ruta_genoma
    # Fallback a resultados/tablas/
    ruta_legacy = os.path.join(RUTA_RESULTADOS_TABLAS, nombre_archivo)
    if os.path.exists(ruta_legacy):
        return ruta_legacy
    return ruta_legacy  # retornar la ruta aunque no exista, para el mensaje de error

if GENOME_BASENAME:
    ARCHIVO_ANALISIS_CODONES = buscar_archivo(f"analisis_codones_{GENOME_BASENAME}.json")
    ARCHIVO_ANALISIS_GENES = buscar_archivo(f"analisis_genes_{GENOME_BASENAME}.json")
else:
    ARCHIVO_ANALISIS_CODONES = os.path.join(RUTA_RESULTADOS_TABLAS, "analisis_codones_completo.json")
    ARCHIVO_ANALISIS_GENES = os.path.join(RUTA_RESULTADOS_TABLAS, "analisis_genes_completo.json")

# Archivos de comparacion - buscar en carpeta del genoma y en tablas legacy
ARCHIVO_COMPARACION = None
dirs_buscar = []
if GENOME_BASENAME:
    dirs_buscar.append(os.path.join(DIRECTORIO_BASE, "resultados", GENOME_BASENAME, "tablas"))
dirs_buscar.append(RUTA_RESULTADOS_TABLAS)
for dir_buscar in dirs_buscar:
    if os.path.exists(dir_buscar):
        for f in os.listdir(dir_buscar):
            if f.startswith("comparacion_") and f.endswith(".json"):
                ARCHIVO_COMPARACION = os.path.join(dir_buscar, f)
                break
    if ARCHIVO_COMPARACION:
        break
if ARCHIVO_COMPARACION is None:
    ARCHIVO_COMPARACION = os.path.join(RUTA_RESULTADOS_TABLAS, "comparacion_ecoli_vs_salmonella.json")

# Configuracion de estilo (con fallback para diferentes versiones)
try:
    plt.style.use('seaborn-v0_8-whitegrid')
except OSError:
    try:
        plt.style.use('seaborn-whitegrid')
    except OSError:
        plt.style.use('ggplot')
sns.set_palette("husl")

# Colores personalizados
COLORES = {
    'primario': '#2E86AB',
    'secundario': '#A23B72',
    'terciario': '#F18F01',
    'cuaternario': '#C73E1D',
    'verde': '#3A7D44',
    'gris': '#6C757D'
}


# FUNCIONES DE UTILIDAD
# =============================================================================

def crear_directorios():
    """Crea los directorios de figuras si no existen."""
    if not os.path.exists(RUTA_RESULTADOS_FIGURAS):
        os.makedirs(RUTA_RESULTADOS_FIGURAS)
        print(f"[INFO] Directorio creado: {RUTA_RESULTADOS_FIGURAS}")


def cargar_datos_json(ruta_archivo):
    """
    Carga datos desde un archivo JSON.

    Args:
        ruta_archivo: Ruta al archivo JSON

    Returns:
        dict: Datos cargados
    """
    if not os.path.exists(ruta_archivo):
        raise FileNotFoundError(f"No se encontro el archivo: {ruta_archivo}")

    with open(ruta_archivo, 'r', encoding='utf-8') as archivo:
        return json.load(archivo)


def guardar_figura(figura, nombre_archivo):
    """
    Guarda una figura en formato PNG.

    Args:
        figura: Objeto figura de matplotlib
        nombre_archivo: Nombre del archivo (sin extension)
    """
    ruta = os.path.join(RUTA_RESULTADOS_FIGURAS, f"{nombre_archivo}.png")
    # DPI reducido a 150 para ahorrar memoria en AWS (suficiente para informes)
    figura.savefig(ruta, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"[OK] Figura guardada: {ruta}")
    plt.close(figura)
    # Limpiar memoria (importante para AWS con poca RAM)
    plt.close('all')
    gc.collect()


# FUNCIONES DE VISUALIZACION - CODONES
# =============================================================================

def graficar_codones_parada(datos_codones):
    """
    Genera graficos de proporciones de codones de parada.

    Args:
        datos_codones: Diccionario con datos de analisis de codones
    """
    print("[INFO] Generando grafico de codones de parada...")

    proporciones = datos_codones['codones_parada']['proporciones_porcentaje']
    conteos = datos_codones['codones_parada']['conteos']
    
    # Usar valores de literatura si están disponibles, sino usar valores observados
    literatura = datos_codones['codones_parada'].get('proporciones_literatura', proporciones)

    codones = list(proporciones.keys())
    valores_obs = [proporciones[c] for c in codones]
    valores_lit = [literatura.get(c, proporciones[c]) for c in codones]
    valores_conteo = [conteos.get(c, 0) for c in codones]

    # Crear figura con dos subgraficos
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Grafico 1: Grafico circular de proporciones observadas
    colores_pastel = [COLORES['primario'], COLORES['secundario'], COLORES['terciario']]

    wedges, texts, autotexts = axes[0].pie(
        valores_obs,
        labels=codones,
        autopct='%1.1f%%',
        colors=colores_pastel,
        explode=(0.02, 0.02, 0.02),
        shadow=True,
        startangle=90
    )

    axes[0].set_title('Proporcion de Codones de Parada\nen el Genoma Completo',
                      fontsize=12, fontweight='bold')

    # Agregar leyenda con conteos
    leyenda_labels = [f'{c}: {conteos[c]:,} codones' for c in codones]
    axes[0].legend(wedges, leyenda_labels, title="Conteos", loc="lower left")

    # Grafico 2: Comparacion con literatura (barras agrupadas)
    x = range(len(codones))
    ancho = 0.35

    barras_obs = axes[1].bar([i - ancho/2 for i in x], valores_obs, ancho,
                              label='Observado (genoma completo)', color=COLORES['primario'])
    barras_lit = axes[1].bar([i + ancho/2 for i in x], valores_lit, ancho,
                              label='Literatura (en genes)', color=COLORES['terciario'])

    axes[1].set_xlabel('Codon de Parada', fontsize=11)
    axes[1].set_ylabel('Proporcion (%)', fontsize=11)
    axes[1].set_title('Comparacion con Literatura Cientifica', fontsize=12, fontweight='bold')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(codones)
    axes[1].legend()
    axes[1].set_ylim(0, 70)

    # Agregar valores sobre las barras
    for barra in barras_obs:
        altura = barra.get_height()
        axes[1].annotate(f'{altura:.1f}%',
                        xy=(barra.get_x() + barra.get_width() / 2, altura),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)

    for barra in barras_lit:
        altura = barra.get_height()
        axes[1].annotate(f'{altura:.1f}%',
                        xy=(barra.get_x() + barra.get_width() / 2, altura),
                        xytext=(0, 3), textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)

    plt.suptitle('Analisis de Codones de Parada - E. coli K-12 MG1655',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    guardar_figura(fig, "codones_parada")


def graficar_codones_inicio(datos_codones):
    """
    Genera grafico de codones de inicio ATG.

    Args:
        datos_codones: Diccionario con datos de analisis de codones
    """
    print("[INFO] Generando grafico de codones de inicio...")

    inicio = datos_codones['codones_inicio']

    fig, ax = plt.subplots(figsize=(10, 6))

    # Datos para el grafico
    categorias = ['Codones ATG\nen genoma', 'Genes\nanotados']
    valores = [inicio['conteo_total'], inicio['genes_anotados_literatura']]
    colores = [COLORES['primario'], COLORES['verde']]

    barras = ax.bar(categorias, valores, color=colores, width=0.5)

    # Agregar valores sobre las barras
    for barra, valor in zip(barras, valores):
        altura = barra.get_height()
        ax.annotate(f'{valor:,}',
                   xy=(barra.get_x() + barra.get_width() / 2, altura),
                   xytext=(0, 5), textcoords="offset points",
                   ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax.set_ylabel('Cantidad', fontsize=11)
    ax.set_title('Codones ATG vs Genes Anotados\nE. coli K-12 MG1655',
                fontsize=14, fontweight='bold')

    # Agregar anotacion explicativa
    ratio = inicio['ratio_atg_por_gen']
    ax.text(0.5, 0.95, f'Ratio: {ratio:.1f} ATG por cada gen anotado',
           transform=ax.transAxes, fontsize=11, ha='center',
           bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax.text(0.5, 0.85, 'No todos los ATG corresponden a inicios de genes funcionales',
           transform=ax.transAxes, fontsize=10, ha='center', style='italic',
           color=COLORES['gris'])

    plt.tight_layout()
    guardar_figura(fig, "codones_inicio_atg")


def graficar_contenido_gc(datos_codones):
    """
    Genera grafico del contenido GC del genoma.

    Args:
        datos_codones: Diccionario con datos de analisis de codones
    """
    print("[INFO] Generando grafico de contenido GC...")

    gc_datos = datos_codones['contenido_gc']
    nucleotidos = gc_datos['conteo_nucleotidos']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Grafico 1: Composicion de nucleotidos
    labels = list(nucleotidos.keys())
    valores = list(nucleotidos.values())
    colores = ['#E74C3C', '#3498DB', '#F39C12', '#27AE60']  # A, T, G, C

    axes[0].pie(valores, labels=labels, autopct='%1.1f%%', colors=colores,
               startangle=90, explode=(0.02, 0.02, 0.02, 0.02))
    axes[0].set_title('Composicion de Nucleotidos', fontsize=12, fontweight='bold')

    # Grafico 2: GC vs AT con comparacion literatura (si está disponible)
    categorias = ['Contenido GC', 'Contenido AT']
    observado = [gc_datos['contenido_gc_porcentaje'], gc_datos['contenido_at_porcentaje']]
    
    # Usar valores de literatura si están disponibles, sino solo mostrar observado
    gc_literatura = gc_datos.get('gc_literatura', None)
    # Asegurar que la variable existe en todos los caminos de ejecución
    barras_lit = []
    
    if gc_literatura is not None:
        literatura = [gc_literatura, 100 - gc_literatura]
        x = range(len(categorias))
        ancho = 0.35

        barras_obs = axes[1].bar([i - ancho/2 for i in x], observado, ancho,
                                 label='Observado', color=COLORES['primario'])
        barras_lit = axes[1].bar([i + ancho/2 for i in x], literatura, ancho,
                                 label='Literatura', color=COLORES['terciario'])

        axes[1].set_ylabel('Porcentaje (%)', fontsize=11)
        axes[1].set_title('Contenido GC vs AT', fontsize=12, fontweight='bold')
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(categorias)
        axes[1].legend()
        axes[1].set_ylim(0, 60)
    else:
        # Si no hay datos de literatura, solo mostrar observado
        x = [0, 1]
        barras_obs = axes[1].bar(x, observado, width=0.6, color=COLORES['primario'])
        axes[1].set_ylabel('Porcentaje (%)', fontsize=11)
        axes[1].set_title('Contenido GC vs AT', fontsize=12, fontweight='bold')
        axes[1].set_xticks(x)
        axes[1].set_xticklabels(categorias)
        axes[1].set_ylim(0, 60)
        
        # Agregar valores en barras
        for barra in barras_obs:
            altura = barra.get_height()
            axes[1].text(barra.get_x() + barra.get_width()/2., altura,
                        f'{altura:.1f}%', ha='center', va='bottom', fontsize=10)

    # Valores sobre barras (solo si hay comparacion con literatura)
    if gc_literatura is not None:
        for barra in list(barras_obs) + list(barras_lit):
            altura = barra.get_height()
            axes[1].annotate(f'{altura:.1f}%',
                            xy=(barra.get_x() + barra.get_width() / 2, altura),
                            xytext=(0, 3), textcoords="offset points",
                            ha='center', va='bottom', fontsize=9)

    plt.suptitle('Contenido de GC - E. coli K-12 MG1655',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    guardar_figura(fig, "contenido_gc")


# FUNCIONES DE VISUALIZACION - GENES
# =============================================================================

def graficar_distribucion_tamanos(datos_genes):
    """
    Genera histograma de distribucion de tamanos de genes.

    Args:
        datos_genes: Diccionario con datos de analisis de genes
    """
    print("[INFO] Generando grafico de distribucion de tamanos...")

    distribucion = datos_genes['distribucion_tamanos']
    estadisticas = datos_genes['estadisticas_generales']['tamano_gen']

    fig, ax = plt.subplots(figsize=(12, 6))

    # Preparar datos - ordenar las claves correctamente
    orden_rangos = ['0-300', '301-600', '601-900', '901-1500', '1501-3000', '>3000']

    # Buscar las claves correctas en el diccionario (pueden variar el formato)
    rangos_ordenados = []
    for patron in orden_rangos:
        for clave in distribucion.keys():
            if patron in clave or clave.startswith(patron):
                rangos_ordenados.append(clave)
                break

    # Si no encontramos todas, usar las claves tal cual vienen
    if len(rangos_ordenados) != len(orden_rangos):
        rangos_ordenados = list(distribucion.keys())

    cantidades = [distribucion[r]['cantidad'] for r in rangos_ordenados]
    porcentajes = [distribucion[r]['porcentaje'] for r in rangos_ordenados]

    # Etiquetas simplificadas para el eje X
    etiquetas = orden_rangos[:len(cantidades)]
    colores = sns.color_palette("Blues_d", len(etiquetas))

    barras = ax.bar(etiquetas, cantidades, color=colores, edgecolor='black', linewidth=0.5)

    # Agregar porcentajes sobre las barras
    for barra, pct in zip(barras, porcentajes):
        altura = barra.get_height()
        ax.annotate(f'{pct:.1f}%',
                   xy=(barra.get_x() + barra.get_width() / 2, altura),
                   xytext=(0, 3), textcoords="offset points",
                   ha='center', va='bottom', fontsize=9)

    ax.set_xlabel('Tamano del Gen (pb)', fontsize=11)
    ax.set_ylabel('Numero de Genes', fontsize=11)
    ax.set_title('Distribucion de Tamanos de Genes\nE. coli K-12 MG1655',
                fontsize=14, fontweight='bold')

    # Agregar lineas de estadisticas
    ax.axhline(y=estadisticas['promedio_pb'], color='red', linestyle='--',
               label=f'Promedio: {estadisticas["promedio_pb"]:.0f} pb', alpha=0.7)

    # Cuadro de estadisticas
    stats_text = (f"Estadisticas de Tamano:\n"
                  f"  Minimo: {estadisticas['minimo_pb']:,} pb\n"
                  f"  Maximo: {estadisticas['maximo_pb']:,} pb\n"
                  f"  Promedio: {estadisticas['promedio_pb']:,.0f} pb\n"
                  f"  Mediana: {estadisticas['mediana_pb']:,.0f} pb")

    ax.text(0.98, 0.95, stats_text, transform=ax.transAxes, fontsize=9,
           verticalalignment='top', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    guardar_figura(fig, "distribucion_tamanos_genes")


def graficar_distribucion_hebras(datos_genes):
    """
    Genera grafico de distribucion de genes por hebra.

    Args:
        datos_genes: Diccionario con datos de analisis de genes
    """
    print("[INFO] Generando grafico de distribucion por hebra...")

    hebras = datos_genes['estadisticas_generales']['distribucion_hebras']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Grafico circular
    valores = [hebras['forward'], hebras['reverse']]
    etiquetas = ['Forward (+)', 'Reverse (-)']
    colores = [COLORES['primario'], COLORES['secundario']]

    axes[0].pie(valores, labels=etiquetas, autopct='%1.1f%%', colors=colores,
               startangle=90, explode=(0.02, 0.02), shadow=True)
    axes[0].set_title('Distribucion de Genes por Hebra', fontsize=12, fontweight='bold')

    # Grafico de barras
    barras = axes[1].bar(etiquetas, valores, color=colores, width=0.5, edgecolor='black')

    for barra, valor in zip(barras, valores):
        altura = barra.get_height()
        axes[1].annotate(f'{valor:,}',
                        xy=(barra.get_x() + barra.get_width() / 2, altura),
                        xytext=(0, 5), textcoords="offset points",
                        ha='center', va='bottom', fontsize=11, fontweight='bold')

    axes[1].set_ylabel('Numero de Genes', fontsize=11)
    axes[1].set_title('Cantidad de Genes por Hebra', fontsize=12, fontweight='bold')

    plt.suptitle('Distribucion de Genes por Hebra - E. coli K-12 MG1655',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    guardar_figura(fig, "distribucion_hebras")


def graficar_comparacion_literatura(datos_genes):
    """
    Genera grafico comparando resultados con literatura.

    Args:
        datos_genes: Diccionario con datos de analisis de genes
    """
    print("[INFO] Generando grafico de comparacion con literatura...")

    comparacion = datos_genes.get('comparacion_literatura', {})
    if not comparacion:
        print("[WARN] No hay datos de comparacion con literatura disponibles")
        return None

    literatura = datos_genes.get('valores_literatura', {})

    fig, ax = plt.subplots(figsize=(10, 6))

    metricas = []
    observados = []
    literatura_vals = []

    # Total CDS
    if 'Total CDS' in comparacion or 'Total genes' in comparacion:
        metricas.append('Total\nCDS')
        total_data = comparacion.get('Total CDS', comparacion.get('Total genes', {}))
        observados.append(total_data.get('observado', 0))
        literatura_vals.append(total_data.get('literatura', 0))

    # Densidad genica
    if 'Densidad genica' in comparacion:
        metricas.append('Densidad\ngenica (%)')
        den_data = comparacion['Densidad genica']
        observados.append(den_data.get('observado', 0))
        literatura_vals.append(den_data.get('literatura', 0))

    # GC en CDS
    if 'GC en CDS' in comparacion:
        metricas.append('GC en\nCDS (%)')
        gc_data = comparacion['GC en CDS']
        observados.append(gc_data.get('observado', 0))
        literatura_vals.append(gc_data.get('literatura', 0))

    # Tamano promedio
    if 'Tamano promedio' in comparacion:
        metricas.append('Tamano\npromedio (pb)')
        tam_data = comparacion['Tamano promedio']
        observados.append(tam_data.get('observado', 0))
        literatura_vals.append(tam_data.get('literatura', 0))

    if not metricas:
        print("[WARN] No hay metricas disponibles para graficar")
        plt.close(fig)
        return None

    # Normalizar para visualizacion (porcentaje del valor de literatura)
    porcentajes_obs = [(o/l)*100 if l != 0 else 100 for o, l in zip(observados, literatura_vals)]

    x = range(len(metricas))
    ancho = 0.35

    # Barras de literatura (100%)
    barras_lit = ax.bar([i - ancho/2 for i in x], [100]*len(metricas), ancho,
                        label='Literatura (referencia)', color=COLORES['gris'], alpha=0.5)

    # Barras observadas
    barras_obs = ax.bar([i + ancho/2 for i in x], porcentajes_obs, ancho,
                        label='Observado', color=COLORES['primario'])

    ax.set_ylabel('Porcentaje respecto a literatura (%)', fontsize=11)
    ax.set_title('Comparacion de Resultados con Literatura Cientifica',
                fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(metricas)
    ax.legend(loc='upper right')
    ax.axhline(y=100, color='red', linestyle='--', alpha=0.5, label='Referencia')
    ax.set_ylim(0, 120)

    # Agregar valores reales como anotacion
    for i, (obs, lit) in enumerate(zip(observados, literatura_vals)):
        ax.annotate(f'Obs: {obs:,.0f}\nLit: {lit:,.0f}',
                   xy=(i, 110), ha='center', fontsize=8,
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7))

    plt.tight_layout()
    guardar_figura(fig, "comparacion_literatura")


def graficar_genes_vs_cds(datos_genes):
    """
    Genera grafico comparando genes totales vs CDS (codificantes).

    Args:
        datos_genes: Diccionario con datos de analisis de genes
    """
    print("[INFO] Generando grafico de genes vs CDS...")

    analisis = datos_genes.get('analisis_genes_vs_cds', {})
    if not analisis:
        print("[WARN] No hay datos de analisis_genes_vs_cds")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # --- Grafico 1: Genes totales vs CDS ---
    ax1 = axes[0]

    categorias = ['Genes totales\nanotados', 'Genes codificantes\n(CDS)', 'Genes no\ncodificantes']
    valores = [
        analisis['genes_totales_anotados'],
        analisis['genes_codificantes_cds'],
        analisis['genes_no_codificantes']
    ]
    colores = [COLORES['gris'], COLORES['primario'], COLORES['secundario']]

    barras = ax1.bar(categorias, valores, color=colores, edgecolor='black', linewidth=0.5)

    for barra, valor in zip(barras, valores):
        altura = barra.get_height()
        ax1.annotate(f'{valor:,}',
                    xy=(barra.get_x() + barra.get_width() / 2, altura),
                    xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax1.set_ylabel('Numero de Genes', fontsize=11)
    ax1.set_title('Genes Totales vs Codificantes (CDS)', fontsize=12, fontweight='bold')

    # --- Grafico 2: Desglose de genes no codificantes ---
    ax2 = axes[1]

    desglose = analisis.get('desglose_rna', {})
    pseudogenes = analisis.get('pseudogenes', 0)

    if desglose:
        etiquetas = list(desglose.keys()) + (['Pseudogenes'] if pseudogenes > 0 else [])
        valores_rna = list(desglose.values()) + ([pseudogenes] if pseudogenes > 0 else [])

        # Descripciones para leyenda
        descripciones = {
            'tRNA': 'ARN transferencia',
            'rRNA': 'ARN ribosomal',
            'ncRNA': 'ARN no codificante',
            'tmRNA': 'ARN tm',
            'Pseudogenes': 'Genes inactivos'
        }

        colores_rna = sns.color_palette("Set2", len(etiquetas))

        wedges, texts, autotexts = ax2.pie(
            valores_rna,
            labels=etiquetas,
            autopct='%1.1f%%',
            colors=colores_rna,
            startangle=90,
            explode=[0.02] * len(etiquetas)
        )

        # Leyenda con descripciones
        leyenda_labels = [f'{e}: {descripciones.get(e, e)}' for e in etiquetas]
        ax2.legend(wedges, leyenda_labels, title="Tipo", loc="lower left", fontsize=8)

    ax2.set_title('Desglose de Genes No Codificantes', fontsize=12, fontweight='bold')

    plt.suptitle('Analisis de Genes Codificantes vs No Codificantes\nE. coli K-12 MG1655',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    guardar_figura(fig, "genes_vs_cds")


def graficar_genes_extremos(datos_genes):
    """
    Genera grafico mostrando los genes extremos (mas largo y mas corto).

    Args:
        datos_genes: Diccionario con datos de analisis de genes
    """
    print("[INFO] Generando grafico de genes extremos...")

    extremos = datos_genes.get('genes_extremos', {})
    if not extremos:
        print("[WARN] No hay datos de genes_extremos")
        return

    gen_largo = extremos['gen_mas_largo']
    gen_corto = extremos['gen_mas_corto']

    fig, ax = plt.subplots(figsize=(12, 7))

    # Datos para las barras
    genes = [f"{gen_corto['locus_tag']}\n({gen_corto['nombre'] or 'sin nombre'})",
             f"{gen_largo['locus_tag']}\n({gen_largo['nombre'] or 'sin nombre'})"]
    tamanos = [gen_corto['longitud_pb'], gen_largo['longitud_pb']]
    colores = [COLORES['secundario'], COLORES['primario']]

    barras = ax.barh(genes, tamanos, color=colores, edgecolor='black', height=0.5)

    # Agregar valores
    for barra, tamano in zip(barras, tamanos):
        ancho = barra.get_width()
        aminoacidos = tamano // 3
        ax.annotate(f'{tamano:,} pb ({aminoacidos:,} aa)',
                   xy=(ancho, barra.get_y() + barra.get_height()/2),
                   xytext=(5, 0), textcoords="offset points",
                   ha='left', va='center', fontsize=10, fontweight='bold')

    ax.set_xlabel('Longitud (pares de bases)', fontsize=11)
    ax.set_title('Genes Extremos: Mas Largo vs Mas Corto\nE. coli K-12 MG1655',
                fontsize=14, fontweight='bold')

    # Cuadro informativo
    info_text = (
        f"GEN MAS LARGO:\n"
        f"  Locus: {gen_largo['locus_tag']}\n"
        f"  Producto: {gen_largo['producto'][:40]}...\n"
        f"  Posicion: {gen_largo['inicio']:,} - {gen_largo['fin']:,}\n\n"
        f"GEN MAS CORTO:\n"
        f"  Locus: {gen_corto['locus_tag']}\n"
        f"  Producto: {gen_corto['producto'][:40]}...\n"
        f"  Posicion: {gen_corto['inicio']:,} - {gen_corto['fin']:,}"
    )

    ax.text(0.98, 0.02, info_text, transform=ax.transAxes, fontsize=8,
           verticalalignment='bottom', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
           family='monospace')

    # Escala logaritmica para mejor visualizacion
    ax.set_xscale('log')
    ax.set_xlim(10, 20000)

    plt.tight_layout()
    guardar_figura(fig, "genes_extremos")


# FUNCIONES DE VISUALIZACION - COMPARACION DE GENOMAS
# =============================================================================

def graficar_comparacion_genomas(datos_comparacion):
    """
    Genera grafico comparando metricas principales entre E. coli y Salmonella.

    Args:
        datos_comparacion: Diccionario con datos de comparacion
    """
    print("[INFO] Generando grafico de comparacion de genomas...")

    metricas = datos_comparacion['metricas_generales']
    ec = metricas['ecoli']
    sal = metricas['salmonella']

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))  # Reducido

    # --- Grafico 1: Tamano del genoma ---
    ax1 = axes[0, 0]
    organismos = ['E. coli K-12', 'Salmonella LT2']
    tamanos = [ec['longitud_genoma_pb'] / 1e6, sal['longitud_genoma_pb'] / 1e6]
    colores = [COLORES['primario'], COLORES['cuaternario']]

    barras = ax1.bar(organismos, tamanos, color=colores, edgecolor='black')
    for barra, tam in zip(barras, tamanos):
        ax1.annotate(f'{tam:.2f} Mb',
                    xy=(barra.get_x() + barra.get_width() / 2, barra.get_height()),
                    xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Megabases (Mb)', fontsize=11)
    ax1.set_title('Tamano del Genoma', fontsize=12, fontweight='bold')

    dif_tamano = sal['longitud_genoma_pb'] - ec['longitud_genoma_pb']
    ax1.text(0.5, 0.02, f'Diferencia: {dif_tamano:,} pb ({dif_tamano/1000:.1f} kb)',
            transform=ax1.transAxes, ha='center', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # --- Grafico 2: Total de genes CDS ---
    ax2 = axes[0, 1]
    genes = [ec['total_genes_cds'], sal['total_genes_cds']]

    barras = ax2.bar(organismos, genes, color=colores, edgecolor='black')
    for barra, gen in zip(barras, genes):
        ax2.annotate(f'{gen:,}',
                    xy=(barra.get_x() + barra.get_width() / 2, barra.get_height()),
                    xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=11, fontweight='bold')
    ax2.set_ylabel('Numero de genes CDS', fontsize=11)
    ax2.set_title('Total de Genes Codificantes', fontsize=12, fontweight='bold')

    dif_genes = sal['total_genes_cds'] - ec['total_genes_cds']
    ax2.text(0.5, 0.02, f'Diferencia: {dif_genes:+,} genes',
            transform=ax2.transAxes, ha='center', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # --- Grafico 3: Contenido GC ---
    ax3 = axes[1, 0]
    gc = [ec['contenido_gc_porcentaje'], sal['contenido_gc_porcentaje']]

    barras = ax3.bar(organismos, gc, color=colores, edgecolor='black')
    for barra, g in zip(barras, gc):
        ax3.annotate(f'{g:.2f}%',
                    xy=(barra.get_x() + barra.get_width() / 2, barra.get_height()),
                    xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=11, fontweight='bold')
    ax3.set_ylabel('Contenido GC (%)', fontsize=11)
    ax3.set_title('Contenido de Guanina + Citosina', fontsize=12, fontweight='bold')
    ax3.set_ylim(0, 60)

    # --- Grafico 4: Densidad genica ---
    ax4 = axes[1, 1]
    densidad = [ec['densidad_genica_porcentaje'], sal['densidad_genica_porcentaje']]

    barras = ax4.bar(organismos, densidad, color=colores, edgecolor='black')
    for barra, d in zip(barras, densidad):
        ax4.annotate(f'{d:.1f}%',
                    xy=(barra.get_x() + barra.get_width() / 2, barra.get_height()),
                    xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Densidad genica (%)', fontsize=11)
    ax4.set_title('Porcentaje del Genoma que Codifica', fontsize=12, fontweight='bold')
    ax4.set_ylim(0, 100)

    plt.suptitle('Comparacion de Genomas: E. coli K-12 vs Salmonella LT2',
                fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    guardar_figura(fig, "comparacion_genomas_metricas")


def graficar_genes_virulencia(datos_comparacion):
    """
    Genera grafico comparando genes de virulencia entre ambos organismos.

    Args:
        datos_comparacion: Diccionario con datos de comparacion
    """
    print("[INFO] Generando grafico de genes de virulencia...")

    virulencia = datos_comparacion['genes_virulencia']
    ec_vir = virulencia['ecoli']
    sal_vir = virulencia['salmonella']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))  # Reducido

    # --- Grafico 1: Total de genes de virulencia ---
    ax1 = axes[0]
    organismos = ['E. coli K-12\n(no patogena)', 'Salmonella LT2\n(patogena)']
    totales = [ec_vir['total'], sal_vir['total']]
    colores = [COLORES['verde'], COLORES['cuaternario']]

    barras = ax1.bar(organismos, totales, color=colores, edgecolor='black', width=0.5)
    for barra, total in zip(barras, totales):
        ax1.annotate(f'{total:,}',
                    xy=(barra.get_x() + barra.get_width() / 2, barra.get_height()),
                    xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=14, fontweight='bold')

    ax1.set_ylabel('Numero de genes', fontsize=11)
    ax1.set_title('Total de Genes Relacionados con Virulencia',
                 fontsize=12, fontweight='bold')

    diferencia = sal_vir['total'] - ec_vir['total']
    ax1.text(0.5, 0.95, f'Salmonella tiene {diferencia:+} genes de virulencia mas',
            transform=ax1.transAxes, ha='center', fontsize=10,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # --- Grafico 2: Categorias de virulencia en Salmonella ---
    ax2 = axes[1]
    categorias = sal_vir['categorias']

    if categorias:
        nombres = list(categorias.keys())
        valores = list(categorias.values())
        colores_cat = sns.color_palette("Reds_r", len(nombres))

        # Ordenar por cantidad
        sorted_data = sorted(zip(valores, nombres, colores_cat), reverse=True)
        valores_ord = [x[0] for x in sorted_data]
        nombres_ord = [x[1] for x in sorted_data]

        barras = ax2.barh(nombres_ord, valores_ord, color=sns.color_palette("Reds_r", len(nombres_ord)))

        for barra, val in zip(barras, valores_ord):
            ax2.annotate(f'{val}',
                        xy=(barra.get_width(), barra.get_y() + barra.get_height()/2),
                        xytext=(3, 0), textcoords="offset points",
                        ha='left', va='center', fontsize=10, fontweight='bold')

        ax2.set_xlabel('Numero de genes', fontsize=11)
        ax2.set_title('Categorias de Genes de Virulencia\nen Salmonella LT2',
                     fontsize=12, fontweight='bold')

    plt.suptitle('Analisis de Genes de Virulencia: E. coli vs Salmonella',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    guardar_figura(fig, "comparacion_virulencia")


def graficar_distancias_intergenicas(datos_comparacion):
    """
    Genera grafico comparando distancias intergenicas (indicador de islas genomicas).

    Args:
        datos_comparacion: Diccionario con datos de comparacion
    """
    print("[INFO] Generando grafico de distancias intergenicas...")

    distancias = datos_comparacion['distancias_intergenicas']
    ec_dist = distancias['ecoli']
    sal_dist = distancias['salmonella']

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))  # Reducido

    # --- Grafico 1: Estadisticas de distancias ---
    ax1 = axes[0]

    metricas = ['Promedio', 'Mediana', 'Maxima']
    ec_vals = [ec_dist['promedio_pb'], ec_dist['mediana_pb'], ec_dist['maxima_pb']]
    sal_vals = [sal_dist['promedio_pb'], sal_dist['mediana_pb'], sal_dist['maxima_pb']]

    x = range(len(metricas))
    ancho = 0.35

    barras_ec = ax1.bar([i - ancho/2 for i in x], ec_vals, ancho,
                        label='E. coli K-12', color=COLORES['primario'])
    barras_sal = ax1.bar([i + ancho/2 for i in x], sal_vals, ancho,
                         label='Salmonella LT2', color=COLORES['cuaternario'])

    ax1.set_ylabel('Distancia (pares de bases)', fontsize=11)
    ax1.set_title('Distancias Intergenicas', fontsize=12, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(metricas)
    ax1.legend()

    # Escala logaritmica para la maxima
    ax1.set_yscale('log')

    # --- Grafico 2: Regiones grandes (posibles islas) ---
    ax2 = axes[1]

    organismos = ['E. coli K-12', 'Salmonella LT2']
    regiones_grandes = [ec_dist['regiones_grandes_5kb'], sal_dist['regiones_grandes_5kb']]
    colores = [COLORES['primario'], COLORES['cuaternario']]

    barras = ax2.bar(organismos, regiones_grandes, color=colores, edgecolor='black', width=0.5)

    for barra, reg in zip(barras, regiones_grandes):
        ax2.annotate(f'{reg}',
                    xy=(barra.get_x() + barra.get_width() / 2, barra.get_height()),
                    xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=14, fontweight='bold')

    ax2.set_ylabel('Numero de regiones', fontsize=11)
    ax2.set_title('Regiones Intergenicas > 5 kb\n(Posibles Islas de Patogenicidad)',
                 fontsize=12, fontweight='bold')

    # Anotacion explicativa
    dif_islas = sal_dist['regiones_grandes_5kb'] - ec_dist['regiones_grandes_5kb']
    if dif_islas > 0:
        texto = f'Salmonella tiene {dif_islas} regiones grandes adicionales\n(pueden contener islas de patogenicidad - SPIs)'
    else:
        texto = 'Distribucion similar de regiones intergenicas'

    ax2.text(0.5, 0.02, texto,
            transform=ax2.transAxes, ha='center', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    plt.suptitle('Distancias Intergenicas: Indicador de Islas Genomicas',
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    guardar_figura(fig, "comparacion_distancias_intergenicas")


def graficar_resumen_comparacion(datos_comparacion):
    """
    Genera un grafico resumen de la comparacion E. coli vs Salmonella.

    Args:
        datos_comparacion: Diccionario con datos de comparacion
    """
    print("[INFO] Generando grafico resumen de comparacion...")

    metricas = datos_comparacion['metricas_generales']
    virulencia = datos_comparacion['genes_virulencia']
    resumen = datos_comparacion.get('resumen_interpretativo', {})

    fig = plt.figure(figsize=(14, 9))  # Reducido para menos memoria

    # Crear grid de subplots
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

    # --- Grafico 1: Comparacion de tamano ---
    ax1 = fig.add_subplot(gs[0, 0])
    organismos = ['E. coli\nK-12', 'Salmonella\nLT2']
    tamanos = [metricas['ecoli']['longitud_genoma_pb']/1e6,
               metricas['salmonella']['longitud_genoma_pb']/1e6]
    colores = [COLORES['primario'], COLORES['cuaternario']]

    ax1.bar(organismos, tamanos, color=colores, edgecolor='black')
    ax1.set_ylabel('Tamano (Mb)')
    ax1.set_title('Tamano del Genoma', fontweight='bold')
    for i, t in enumerate(tamanos):
        ax1.text(i, t + 0.05, f'{t:.2f}', ha='center', fontweight='bold')

    # --- Grafico 2: Comparacion de genes ---
    ax2 = fig.add_subplot(gs[0, 1])
    genes = [metricas['ecoli']['total_genes_cds'],
             metricas['salmonella']['total_genes_cds']]

    ax2.bar(organismos, genes, color=colores, edgecolor='black')
    ax2.set_ylabel('Genes CDS')
    ax2.set_title('Total Genes Codificantes', fontweight='bold')
    for i, g in enumerate(genes):
        ax2.text(i, g + 50, f'{g:,}', ha='center', fontweight='bold')

    # --- Grafico 3: Genes de virulencia ---
    ax3 = fig.add_subplot(gs[0, 2])
    vir = [virulencia['ecoli']['total'], virulencia['salmonella']['total']]

    ax3.bar(organismos, vir, color=[COLORES['verde'], COLORES['cuaternario']], edgecolor='black')
    ax3.set_ylabel('Genes de virulencia')
    ax3.set_title('Factores de Patogenicidad', fontweight='bold')
    for i, v in enumerate(vir):
        ax3.text(i, v + 2, f'{v}', ha='center', fontweight='bold')

    # --- Grafico 4: Cuadro comparativo (tabla simple) ---
    ax4 = fig.add_subplot(gs[1, :])
    ax4.axis('off')

    # Crear tabla de diferencias usando una tabla real de matplotlib
    diferencias = metricas['diferencias']

    # Datos para la tabla
    tabla_data = [
        ['TAMANO DEL GENOMA', f'Salmonella tiene {diferencias["longitud_pb"]:,} pb adicionales ({diferencias["longitud_porcentaje"]:.1f}% mas)'],
        ['GENES CODIFICANTES', f'Salmonella tiene +{diferencias["total_genes"]} genes CDS'],
        ['GENES DE VIRULENCIA', f'Salmonella tiene +{virulencia["diferencia_total"]} genes de patogenicidad'],
    ]

    tabla = ax4.table(cellText=tabla_data,
                      colLabels=['Metrica', 'Diferencia'],
                      loc='upper center',
                      cellLoc='left',
                      colWidths=[0.3, 0.6])

    tabla.auto_set_font_size(False)
    tabla.set_fontsize(9)
    tabla.scale(1, 1.8)

    # Estilo de la tabla
    for i in range(len(tabla_data) + 1):
        for j in range(2):
            cell = tabla[(i, j)]
            if i == 0:  # Header
                cell.set_facecolor('#2E86AB')
                cell.set_text_props(color='white', fontweight='bold')
            else:
                cell.set_facecolor('#f8f9fa' if i % 2 == 0 else 'white')

    # Texto explicativo debajo
    explicacion = (
        "POR QUE SALMONELLA ES PATOGENA:\n"
        "- Posee sistemas de secrecion tipo III (T3SS)\n"
        "- Tiene islas de patogenicidad (SPIs)\n"
        "- Puede sobrevivir dentro de macrofagos\n\n"
        "POR QUE E. coli K-12 NO ES PATOGENA:\n"
        "- Cepa de laboratorio domesticada desde 1922\n"
        "- Ha perdido genes de virulencia\n"
        "- Es el organismo modelo mas seguro"
    )

    ax4.text(0.5, 0.15, explicacion, transform=ax4.transAxes,
            fontsize=8, verticalalignment='center', horizontalalignment='center',
            bbox=dict(boxstyle='round', facecolor='lightyellow', edgecolor='orange', alpha=0.8))

    plt.suptitle('Resumen Comparativo: E. coli K-12 vs Salmonella LT2\n'
                 'Entendiendo las diferencias entre comensal y patogeno',
                fontsize=14, fontweight='bold', y=0.98)

    guardar_figura(fig, "resumen_comparacion_ecoli_salmonella")


def graficar_resumen_general(datos_codones, datos_genes):
    """
    Genera un grafico resumen con las metricas principales.

    Args:
        datos_codones: Diccionario con datos de codones
        datos_genes: Diccionario con datos de genes
    """
    print("[INFO] Generando grafico resumen general...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))  # Reducido para menos memoria

    # --- Subgrafico 1: Metricas principales ---
    ax1 = axes[0, 0]

    metricas = ['Longitud\ngenoma (Mpb)', 'Total\ngenes', 'Contenido\nGC (%)', 'Densidad\ngenica (%)']
    valores = [
        datos_codones['codones_inicio']['longitud_genoma_pb'] / 1e6,
        datos_genes['estadisticas_generales']['total_genes'],
        datos_codones['contenido_gc']['contenido_gc_porcentaje'],
        datos_genes['estadisticas_generales']['densidad_genica_porcentaje']
    ]

    colores_metricas = [COLORES['primario'], COLORES['verde'],
                        COLORES['terciario'], COLORES['secundario']]

    barras = ax1.bar(metricas, valores, color=colores_metricas)
    ax1.set_title('Metricas Principales del Genoma', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Valor', fontsize=10)

    for barra, valor in zip(barras, valores):
        altura = barra.get_height()
        formato = f'{valor:.2f}' if valor < 100 else f'{valor:,.0f}'
        ax1.annotate(formato, xy=(barra.get_x() + barra.get_width() / 2, altura),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    # --- Subgrafico 2: Codones de parada ---
    ax2 = axes[0, 1]

    proporciones = datos_codones['codones_parada']['proporciones_porcentaje']
    codones = list(proporciones.keys())
    valores_parada = [proporciones[c] for c in codones]

    colores_codones = [COLORES['primario'], COLORES['secundario'], COLORES['terciario']]
    ax2.pie(valores_parada, labels=codones, autopct='%1.1f%%', colors=colores_codones,
           startangle=90, explode=(0.02, 0.02, 0.02))
    ax2.set_title('Codones de Parada', fontsize=12, fontweight='bold')

    # --- Subgrafico 3: Distribucion de tamanos ---
    ax3 = axes[1, 0]

    distribucion = datos_genes['distribucion_tamanos']
    etiquetas_dist = ['0-300', '301-600', '601-900', '901-1500', '1501-3000', '>3000']

    # Ordenar claves correctamente
    rangos_ordenados = []
    for patron in etiquetas_dist:
        for clave in distribucion.keys():
            if patron in clave or clave.startswith(patron):
                rangos_ordenados.append(clave)
                break

    if len(rangos_ordenados) != len(etiquetas_dist):
        rangos_ordenados = list(distribucion.keys())

    cantidades = [distribucion[r]['cantidad'] for r in rangos_ordenados]

    ax3.bar(etiquetas_dist[:len(cantidades)], cantidades, color=sns.color_palette("Blues_d", len(cantidades)))
    ax3.set_xlabel('Tamano (pb)', fontsize=10)
    ax3.set_ylabel('Genes', fontsize=10)
    ax3.set_title('Distribucion de Tamanos', fontsize=12, fontweight='bold')
    ax3.tick_params(axis='x', rotation=45)

    # --- Subgrafico 4: Hebras ---
    ax4 = axes[1, 1]

    hebras = datos_genes['estadisticas_generales']['distribucion_hebras']
    ax4.pie([hebras['forward'], hebras['reverse']],
           labels=['Forward (+)', 'Reverse (-)'],
           autopct='%1.1f%%', colors=[COLORES['primario'], COLORES['secundario']],
           startangle=90)
    ax4.set_title('Genes por Hebra', fontsize=12, fontweight='bold')

    plt.suptitle('Resumen del Analisis del Genoma de E. coli K-12 MG1655',
                fontsize=16, fontweight='bold', y=1.02)
    plt.tight_layout()

    guardar_figura(fig, "resumen_general")


# FUNCION PRINCIPAL
# =============================================================================

def main():
    """Funcion principal que genera todas las visualizaciones."""

    print("\n" + "=" * 70)
    print("GENERACION DE VISUALIZACIONES - GENOMAS BACTERIANOS")
    print("=" * 70)
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Crear directorios
    crear_directorios()

    # Cargar datos
    print("\n[INFO] Cargando datos de analisis...")

    datos_codones = None
    datos_genes = None
    datos_comparacion = None

    try:
        datos_codones = cargar_datos_json(ARCHIVO_ANALISIS_CODONES)
        print(f"[OK] Datos de codones cargados")
    except FileNotFoundError:
        print(f"[WARN] No se encontro: {ARCHIVO_ANALISIS_CODONES}")
        print("       Ejecuta analisis_codones.py para generar graficos de codones")

    try:
        datos_genes = cargar_datos_json(ARCHIVO_ANALISIS_GENES)
        print(f"[OK] Datos de genes cargados")
    except FileNotFoundError:
        print(f"[WARN] No se encontro: {ARCHIVO_ANALISIS_GENES}")
        print("       Ejecuta analisis_genes.py para generar graficos de genes")

    # Cargar datos de comparacion (si existen)
    try:
        datos_comparacion = cargar_datos_json(ARCHIVO_COMPARACION)
        print(f"[OK] Datos de comparacion E. coli vs Salmonella cargados")
    except FileNotFoundError:
        print(f"[INFO] No se encontro archivo de comparacion")
        print("       Ejecuta comparar_genomas.py para generar graficos comparativos")

    # Verificar que hay al menos algunos datos
    if datos_codones is None and datos_genes is None and datos_comparacion is None:
        print("\n[ERROR] No se encontraron datos para generar visualizaciones")
        print("        Ejecuta primero los scripts de analisis:")
        print("        - python analisis_codones.py")
        print("        - python analisis_genes.py")
        print("        - python comparar_genomas.py")
        return

    # Generar visualizaciones
    print("\n" + "=" * 70)
    print("GENERANDO GRAFICOS")
    print("=" * 70)

    figuras_generadas = []

    # Funcion auxiliar para ejecutar graficos con manejo de errores
    def ejecutar_grafico(funcion, nombre_archivo, *args):
        """Ejecuta una funcion de grafico con manejo de errores."""
        try:
            funcion(*args)
            figuras_generadas.append(nombre_archivo)
            return True
        except Exception as e:
            print(f"[ERROR] Fallo al generar {nombre_archivo}: {e}")
            errores.append((nombre_archivo, str(e)))
            plt.close('all')
            gc.collect()
            return False

    errores = []

    # Graficos de codones (si hay datos)
    if datos_codones:
        print("\n[SECCION] Graficos de analisis de codones...")
        ejecutar_grafico(graficar_codones_parada, "codones_parada.png", datos_codones)
        ejecutar_grafico(graficar_codones_inicio, "codones_inicio_atg.png", datos_codones)
        ejecutar_grafico(graficar_contenido_gc, "contenido_gc.png", datos_codones)

    # Graficos de genes (si hay datos)
    if datos_genes:
        print("\n[SECCION] Graficos de analisis de genes...")
        ejecutar_grafico(graficar_distribucion_tamanos, "distribucion_tamanos_genes.png", datos_genes)
        ejecutar_grafico(graficar_distribucion_hebras, "distribucion_hebras.png", datos_genes)
        ejecutar_grafico(graficar_comparacion_literatura, "comparacion_literatura.png", datos_genes)
        ejecutar_grafico(graficar_genes_vs_cds, "genes_vs_cds.png", datos_genes)
        ejecutar_grafico(graficar_genes_extremos, "genes_extremos.png", datos_genes)

    # Grafico resumen (si hay ambos datos)
    if datos_codones and datos_genes:
        print("\n[SECCION] Grafico resumen general...")
        ejecutar_grafico(graficar_resumen_general, "resumen_general.png", datos_codones, datos_genes)

    # Graficos de comparacion E. coli vs Salmonella (si hay datos)
    if datos_comparacion:
        print("\n[SECCION] Graficos de comparacion E. coli vs Salmonella...")
        ejecutar_grafico(graficar_comparacion_genomas, "comparacion_genomas_metricas.png", datos_comparacion)
        ejecutar_grafico(graficar_genes_virulencia, "comparacion_virulencia.png", datos_comparacion)
        ejecutar_grafico(graficar_distancias_intergenicas, "comparacion_distancias_intergenicas.png", datos_comparacion)
        ejecutar_grafico(graficar_resumen_comparacion, "resumen_comparacion_ecoli_salmonella.png", datos_comparacion)

    # Resumen
    print("\n" + "=" * 70)
    print("FIGURAS GENERADAS")
    print("=" * 70)
    print(f"  Directorio: {RUTA_RESULTADOS_FIGURAS}")
    print(f"  Total de figuras: {len(figuras_generadas)}")
    print("\n  Archivos generados:")

    # Mostrar figuras por categoria
    if datos_codones:
        print("\n  [Analisis de Codones]")
        print("    - codones_parada.png          (proporciones de codones de parada)")
        print("    - codones_inicio_atg.png      (codones ATG vs genes)")
        print("    - contenido_gc.png            (composicion de nucleotidos)")

    if datos_genes:
        print("\n  [Analisis de Genes]")
        print("    - distribucion_tamanos_genes.png (histograma de tamanos)")
        print("    - distribucion_hebras.png     (genes por hebra + / -)")
        print("    - comparacion_literatura.png  (validacion con fuentes)")
        print("    - genes_vs_cds.png            (codificantes vs no codificantes)")
        print("    - genes_extremos.png          (gen mas largo y mas corto)")

    if datos_codones and datos_genes:
        print("\n  [Resumen]")
        print("    - resumen_general.png         (dashboard completo)")

    if datos_comparacion:
        print("\n  [Comparacion E. coli vs Salmonella]")
        print("    - comparacion_genomas_metricas.png     (metricas lado a lado)")
        print("    - comparacion_virulencia.png          (genes de patogenicidad)")
        print("    - comparacion_distancias_intergenicas.png (islas genomicas)")
        print("    - resumen_comparacion_ecoli_salmonella.png (resumen comparativo)")

    # Mostrar errores si los hubo
    if errores:
        print("\n  [ERRORES]")
        for archivo, error in errores:
            print(f"    - {archivo}: {error}")

    print("\n" + "=" * 70)
    if errores:
        print(f"[WARN] {len(figuras_generadas)} visualizaciones generadas, {len(errores)} con errores")
    else:
        print(f"[OK] {len(figuras_generadas)} visualizaciones generadas exitosamente")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()

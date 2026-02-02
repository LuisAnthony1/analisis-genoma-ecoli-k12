#!/usr/bin/env python3
"""
Visualizaciones del Analisis del Genoma de E. coli K-12 MG1655

Este script genera graficos profesionales para ilustrar los resultados
del analisis del genoma de E. coli K-12 MG1655.

Graficos generados:
- Distribucion de tamanos de genes (histograma)
- Proporciones de codones de parada (grafico circular y barras)
- Contenido GC por gen (histograma)
- Distribucion de genes por hebra (barras)
- Comparacion de resultados con literatura cientifica

Fecha: 2026
Proyecto: Analisis del genoma de E. coli K-12 MG1655
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import json
import os
from datetime import datetime


# CONFIGURACION
# =============================================================================

# Rutas de archivos
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUTA_RESULTADOS_TABLAS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")
RUTA_RESULTADOS_FIGURAS = os.path.join(DIRECTORIO_BASE, "resultados", "figuras")

# Archivos de entrada (generados por los scripts de analisis)
ARCHIVO_ANALISIS_CODONES = os.path.join(RUTA_RESULTADOS_TABLAS, "analisis_codones_completo.json")
ARCHIVO_ANALISIS_GENES = os.path.join(RUTA_RESULTADOS_TABLAS, "analisis_genes_completo.json")

# Configuracion de estilo
plt.style.use('seaborn-v0_8-whitegrid')
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
    figura.savefig(ruta, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"[OK] Figura guardada: {ruta}")
    plt.close(figura)


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
    literatura = datos_codones['codones_parada']['proporciones_literatura']

    codones = list(proporciones.keys())
    valores_obs = [proporciones[c] for c in codones]
    valores_lit = [literatura[c] for c in codones]
    valores_conteo = [conteos[c] for c in codones]

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

    # Grafico 2: GC vs AT con comparacion literatura
    categorias = ['Contenido GC', 'Contenido AT']
    observado = [gc_datos['contenido_gc_porcentaje'], gc_datos['contenido_at_porcentaje']]
    literatura = [gc_datos['gc_literatura'], 100 - gc_datos['gc_literatura']]

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

    # Valores sobre barras
    for barra in barras_obs + barras_lit:
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

    # Preparar datos
    rangos = list(distribucion.keys())
    cantidades = [distribucion[r]['cantidad'] for r in rangos]
    porcentajes = [distribucion[r]['porcentaje'] for r in rangos]

    # Simplificar etiquetas
    etiquetas = ['0-300', '301-600', '601-900', '901-1500', '1501-3000', '>3000']
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

    comparacion = datos_genes['comparacion_literatura']
    literatura = datos_genes['valores_literatura']

    fig, ax = plt.subplots(figsize=(10, 6))

    metricas = ['Total\ngenes', 'Densidad\ngenica (%)', 'GC en\nCDS (%)', 'Tamano\npromedio (pb)']

    observados = [
        comparacion['Total genes']['observado'],
        comparacion['Densidad genica']['observado'],
        comparacion['GC en CDS']['observado'],
        comparacion['Tamano promedio']['observado']
    ]

    literatura_vals = [
        comparacion['Total genes']['literatura'],
        comparacion['Densidad genica']['literatura'],
        comparacion['GC en CDS']['literatura'],
        comparacion['Tamano promedio']['literatura']
    ]

    # Normalizar para visualizacion (porcentaje del valor de literatura)
    porcentajes_obs = [(o/l)*100 if l != 0 else 100 for o, l in zip(observados, literatura_vals)]

    x = range(len(metricas))
    ancho = 0.35

    # Barras de literatura (100%)
    barras_lit = ax.bar([i - ancho/2 for i in x], [100]*4, ancho,
                        label='Literatura (referencia)', color=COLORES['gris'], alpha=0.5)

    # Barras observadas
    barras_obs = ax.bar([i + ancho/2 for i in x], porcentajes_obs, ancho,
                        label='Observado', color=COLORES['primario'])

    ax.set_ylabel('Porcentaje respecto a literatura (%)', fontsize=11)
    ax.set_title('Comparacion de Resultados con Literatura Cientifica\nE. coli K-12 MG1655',
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


def graficar_resumen_general(datos_codones, datos_genes):
    """
    Genera un grafico resumen con las metricas principales.

    Args:
        datos_codones: Diccionario con datos de codones
        datos_genes: Diccionario con datos de genes
    """
    print("[INFO] Generando grafico resumen general...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

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
    rangos_keys = list(distribucion.keys())
    cantidades = [distribucion[r]['cantidad'] for r in rangos_keys]

    ax3.bar(etiquetas_dist, cantidades, color=sns.color_palette("Blues_d", 6))
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

    print("\n" + "=" * 60)
    print("GENERACION DE VISUALIZACIONES - E. coli K-12 MG1655")
    print("=" * 60)
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Crear directorios
    crear_directorios()

    # Cargar datos
    print("\n[INFO] Cargando datos de analisis...")

    try:
        datos_codones = cargar_datos_json(ARCHIVO_ANALISIS_CODONES)
        print(f"[OK] Datos de codones cargados")
    except FileNotFoundError:
        print(f"[ERROR] No se encontro: {ARCHIVO_ANALISIS_CODONES}")
        print("[INFO] Ejecuta primero analisis_codones.py")
        return

    try:
        datos_genes = cargar_datos_json(ARCHIVO_ANALISIS_GENES)
        print(f"[OK] Datos de genes cargados")
    except FileNotFoundError:
        print(f"[ERROR] No se encontro: {ARCHIVO_ANALISIS_GENES}")
        print("[INFO] Ejecuta primero analisis_genes.py")
        return

    # Generar visualizaciones
    print("\n" + "=" * 60)
    print("GENERANDO GRAFICOS")
    print("=" * 60)

    # Graficos de codones
    graficar_codones_parada(datos_codones)
    graficar_codones_inicio(datos_codones)
    graficar_contenido_gc(datos_codones)

    # Graficos de genes
    graficar_distribucion_tamanos(datos_genes)
    graficar_distribucion_hebras(datos_genes)
    graficar_comparacion_literatura(datos_genes)

    # Grafico resumen
    graficar_resumen_general(datos_codones, datos_genes)

    # Resumen
    print("\n" + "=" * 60)
    print("FIGURAS GENERADAS")
    print("=" * 60)
    print(f"  Directorio: {RUTA_RESULTADOS_FIGURAS}")
    print("  Archivos:")
    print("    - codones_parada.png")
    print("    - codones_inicio_atg.png")
    print("    - contenido_gc.png")
    print("    - distribucion_tamanos_genes.png")
    print("    - distribucion_hebras.png")
    print("    - comparacion_literatura.png")
    print("    - resumen_general.png")
    print("=" * 60)

    print("\n[OK] Todas las visualizaciones generadas exitosamente\n")


if __name__ == "__main__":
    main()

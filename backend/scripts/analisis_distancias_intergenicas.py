#!/usr/bin/env python3
"""
Análisis de Distancias Intergénicas del Genoma de E. coli K-12 MG1655

Este script calcula las distancias entre genes consecutivos en el genoma
de E. coli K-12 MG1655, analizando las regiones intergénicas.

Funcionalidades:
- Calcular distancias entre genes consecutivos
- Identificar genes solapados (distancia negativa)
- Analizar estadísticas de espaciamiento génico
- Detectar regiones intergénicas grandes (posibles promotores, operones)
- Comparar distancias en hebra forward (+) vs reverse (-)
- Exportar resultados en formato CSV y JS

Fecha: 2026
Proyecto: Análisis del genoma de E. coli K-12 MG1655
"""

from Bio import SeqIO
import os
import json
import csv
from datetime import datetime
import statistics
from collections import Counter


# CONFIGURACIÓN
# =============================================================================

# Rutas de archivos
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_BASE, "datos", "crudo")
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")
ARCHIVO_GENBANK = os.path.join(RUTA_DATOS_CRUDO, "ecoli_k12.gb")

# Crear directorio de resultados si no existe
os.makedirs(RUTA_RESULTADOS, exist_ok=True)

# Umbrales de análisis
UMBRAL_SOLAPAMIENTO = 0  # Genes con distancia < 0 están solapados
UMBRAL_REGION_GRANDE = 500  # Regiones intergénicas > 500 pb son "grandes"


# FUNCIONES DE EXTRACCIÓN
# =============================================================================

def extraer_genes_ordenados(archivo_genbank):
    """
    Extrae todos los genes CDS del GenBank y los ordena por posición.
    
    Args:
        archivo_genbank: Ruta al archivo GenBank
    
    Returns:
        tuple: (lista de genes ordenados, longitud del genoma)
    """
    print("\n" + "=" * 60)
    print("EXTRAYENDO GENES DEL GENOMA")
    print("=" * 60)
    
    # Leer el genoma
    registro = SeqIO.read(archivo_genbank, "genbank")
    longitud_genoma = len(registro.seq)
    
    print(f"Archivo: {os.path.basename(archivo_genbank)}")
    print(f"Organismo: {registro.description}")
    print(f"Longitud del genoma: {longitud_genoma:,} pb")
    
    genes = []
    
    # Extraer información de cada gen CDS
    for feature in registro.features:
        if feature.type == "CDS":
            inicio = int(feature.location.start)
            fin = int(feature.location.end)
            hebra = feature.location.strand  # 1 = forward, -1 = reverse
            
            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
            nombre_gen = feature.qualifiers.get("gene", [""])[0]
            producto = feature.qualifiers.get("product", [""])[0]
            
            gen_info = {
                "locus_tag": locus_tag,
                "nombre_gen": nombre_gen,
                "inicio": inicio,
                "fin": fin,
                "longitud": fin - inicio,
                "hebra": "+" if hebra == 1 else "-",
                "producto": producto
            }
            
            genes.append(gen_info)
    
    # Ordenar genes por posición de inicio
    genes_ordenados = sorted(genes, key=lambda x: x["inicio"])
    
    print(f"[OK] Genes extraídos y ordenados: {len(genes_ordenados):,} CDS")
    
    return genes_ordenados, longitud_genoma


# FUNCIONES DE ANÁLISIS
# =============================================================================

def calcular_distancias_intergenicas(genes):
    """
    Calcula la distancia entre genes consecutivos.
    
    La distancia intergénica se define como:
    - Positiva: espacio entre el fin del gen N y el inicio del gen N+1
    - Cero: genes adyacentes sin espaciamiento
    - Negativa: genes solapados
    
    Args:
        genes: Lista de genes ordenados por posición
    
    Returns:
        list: Lista de distancias intergénicas con información adicional
    """
    print("\n" + "=" * 60)
    print("CALCULANDO DISTANCIAS INTERGÉNICAS")
    print("=" * 60)
    
    distancias = []
    
    for i in range(len(genes) - 1):
        gen_actual = genes[i]
        gen_siguiente = genes[i + 1]
        
        # Distancia = inicio del siguiente - fin del actual
        distancia = gen_siguiente["inicio"] - gen_actual["fin"]
        
        distancia_info = {
            "gen1_locus": gen_actual["locus_tag"],
            "gen1_nombre": gen_actual["nombre_gen"],
            "gen1_fin": gen_actual["fin"],
            "gen1_hebra": gen_actual["hebra"],
            "gen2_locus": gen_siguiente["locus_tag"],
            "gen2_nombre": gen_siguiente["nombre_gen"],
            "gen2_inicio": gen_siguiente["inicio"],
            "gen2_hebra": gen_siguiente["hebra"],
            "distancia_pb": distancia,
            "tipo": clasificar_distancia(distancia),
            "misma_hebra": gen_actual["hebra"] == gen_siguiente["hebra"]
        }
        
        distancias.append(distancia_info)
    
    print(f"[OK] Distancias calculadas: {len(distancias):,} pares de genes")
    
    return distancias


def clasificar_distancia(distancia):
    """
    Clasifica el tipo de región intergénica según la distancia.
    
    Args:
        distancia: Distancia en pares de bases
    
    Returns:
        str: Clasificación de la distancia
    """
    if distancia < 0:
        return "Solapado"
    elif distancia == 0:
        return "Adyacente"
    elif distancia < 50:
        return "Muy cercano (< 50 pb)"
    elif distancia < 100:
        return "Cercano (50-100 pb)"
    elif distancia < 200:
        return "Moderado (100-200 pb)"
    elif distancia < 500:
        return "Normal (200-500 pb)"
    else:
        return "Grande (> 500 pb)"


def analizar_estadisticas_distancias(distancias):
    """
    Calcula estadísticas generales de las distancias intergénicas.
    
    Args:
        distancias: Lista de distancias intergénicas
    
    Returns:
        dict: Estadísticas completas
    """
    print("\n" + "=" * 60)
    print("ESTADÍSTICAS DE DISTANCIAS INTERGÉNICAS")
    print("=" * 60)
    
    # Extraer solo los valores numéricos
    valores_distancias = [d["distancia_pb"] for d in distancias]
    
    # Separar por categorías
    solapados = [d for d in valores_distancias if d < 0]
    adyacentes = [d for d in valores_distancias if d == 0]
    espaciados = [d for d in valores_distancias if d > 0]
    grandes = [d for d in valores_distancias if d >= UMBRAL_REGION_GRANDE]
    
    # Estadísticas generales
    total_pares = len(valores_distancias)
    
    estadisticas = {
        "total_pares_genes": total_pares,
        "genes_solapados": len(solapados),
        "genes_adyacentes": len(adyacentes),
        "genes_espaciados": len(espaciados),
        "regiones_grandes": len(grandes),
        "porcentaje_solapados": round((len(solapados) / total_pares) * 100, 2),
        "porcentaje_adyacentes": round((len(adyacentes) / total_pares) * 100, 2),
        "porcentaje_espaciados": round((len(espaciados) / total_pares) * 100, 2),
        "porcentaje_grandes": round((len(grandes) / total_pares) * 100, 2),
    }
    
    # Estadísticas solo de regiones espaciadas (distancia > 0)
    if espaciados:
        estadisticas.update({
            "distancia_minima": min(espaciados),
            "distancia_maxima": max(espaciados),
            "distancia_promedio": round(statistics.mean(espaciados), 2),
            "distancia_mediana": statistics.median(espaciados),
            "desviacion_std": round(statistics.stdev(espaciados), 2),
        })
    
    # Estadísticas de solapamiento
    if solapados:
        estadisticas["solapamiento_maximo"] = min(solapados)  # El más negativo
        estadisticas["solapamiento_promedio"] = round(statistics.mean(solapados), 2)
    
    # Imprimir resultados
    print(f"\nTotal de pares analizados: {total_pares:,}")
    print(f"\nCLASIFICACIÓN:")
    print(f"  • Genes solapados:      {len(solapados):5,} ({estadisticas['porcentaje_solapados']:5.2f}%)")
    print(f"  • Genes adyacentes:     {len(adyacentes):5,} ({estadisticas['porcentaje_adyacentes']:5.2f}%)")
    print(f"  • Genes espaciados:     {len(espaciados):5,} ({estadisticas['porcentaje_espaciados']:5.2f}%)")
    print(f"  • Regiones grandes:     {len(grandes):5,} ({estadisticas['porcentaje_grandes']:5.2f}%)")
    
    if espaciados:
        print(f"\nESTADÍSTICAS DE ESPACIAMIENTO (distancias > 0):")
        print(f"  • Distancia mínima:     {estadisticas['distancia_minima']:,} pb")
        print(f"  • Distancia máxima:     {estadisticas['distancia_maxima']:,} pb")
        print(f"  • Distancia promedio:   {estadisticas['distancia_promedio']:,.2f} pb")
        print(f"  • Distancia mediana:    {estadisticas['distancia_mediana']:,} pb")
        print(f"  • Desviación estándar:  {estadisticas['desviacion_std']:,.2f} pb")
    
    if solapados:
        print(f"\nESTADÍSTICAS DE SOLAPAMIENTO:")
        print(f"  • Solapamiento máximo:  {abs(estadisticas['solapamiento_maximo']):,} pb")
        print(f"  • Solapamiento promedio: {abs(estadisticas['solapamiento_promedio']):,.2f} pb")
    
    return estadisticas


def analizar_por_hebra(distancias):
    """
    Analiza las distancias según si los genes están en la misma hebra o no.
    
    Args:
        distancias: Lista de distancias intergénicas
    
    Returns:
        dict: Estadísticas por hebra
    """
    print("\n" + "=" * 60)
    print("ANÁLISIS POR HEBRA")
    print("=" * 60)
    
    misma_hebra = [d["distancia_pb"] for d in distancias if d["misma_hebra"]]
    diferente_hebra = [d["distancia_pb"] for d in distancias if not d["misma_hebra"]]
    
    # Filtrar solo distancias positivas para promedios
    misma_hebra_pos = [d for d in misma_hebra if d > 0]
    diferente_hebra_pos = [d for d in diferente_hebra if d > 0]
    
    estadisticas_hebra = {
        "misma_hebra_total": len(misma_hebra),
        "diferente_hebra_total": len(diferente_hebra),
    }
    
    if misma_hebra_pos:
        estadisticas_hebra["misma_hebra_promedio"] = round(statistics.mean(misma_hebra_pos), 2)
        estadisticas_hebra["misma_hebra_mediana"] = statistics.median(misma_hebra_pos)
    
    if diferente_hebra_pos:
        estadisticas_hebra["diferente_hebra_promedio"] = round(statistics.mean(diferente_hebra_pos), 2)
        estadisticas_hebra["diferente_hebra_mediana"] = statistics.median(diferente_hebra_pos)
    
    print(f"\nGenes en la MISMA hebra:")
    print(f"  • Total de pares:       {len(misma_hebra):,}")
    if misma_hebra_pos:
        print(f"  • Distancia promedio:   {estadisticas_hebra['misma_hebra_promedio']:,.2f} pb")
        print(f"  • Distancia mediana:    {estadisticas_hebra['misma_hebra_mediana']:,} pb")
    
    print(f"\nGenes en DIFERENTE hebra:")
    print(f"  • Total de pares:       {len(diferente_hebra):,}")
    if diferente_hebra_pos:
        print(f"  • Distancia promedio:   {estadisticas_hebra['diferente_hebra_promedio']:,.2f} pb")
        print(f"  • Distancia mediana:    {estadisticas_hebra['diferente_hebra_mediana']:,} pb")
    
    return estadisticas_hebra


def encontrar_regiones_grandes(distancias, umbral=UMBRAL_REGION_GRANDE):
    """
    Identifica las regiones intergénicas más grandes.
    
    Args:
        distancias: Lista de distancias intergénicas
        umbral: Umbral mínimo para considerar una región como "grande"
    
    Returns:
        list: Regiones grandes ordenadas por tamaño
    """
    print("\n" + "=" * 60)
    print(f"REGIONES INTERGÉNICAS GRANDES (> {umbral} pb)")
    print("=" * 60)
    
    regiones_grandes = [d for d in distancias if d["distancia_pb"] >= umbral]
    regiones_grandes_ordenadas = sorted(regiones_grandes, 
                                       key=lambda x: x["distancia_pb"], 
                                       reverse=True)
    
    print(f"\nTotal de regiones grandes: {len(regiones_grandes):,}")
    print(f"\nTop 10 regiones más grandes:")
    print(f"{'#':<4} {'Distancia':<12} {'Gen 1':<15} {'Gen 2':<15} {'Hebras':<10}")
    print("-" * 60)
    
    for i, region in enumerate(regiones_grandes_ordenadas[:10], 1):
        gen1 = region['gen1_nombre'] or region['gen1_locus']
        gen2 = region['gen2_nombre'] or region['gen2_locus']
        hebras = f"{region['gen1_hebra']} → {region['gen2_hebra']}"
        
        print(f"{i:<4} {region['distancia_pb']:>10,} pb  "
              f"{gen1:<15} {gen2:<15} {hebras:<10}")
    
    return regiones_grandes_ordenadas


def analizar_distribucion_tipos(distancias):
    """
    Analiza la distribución de tipos de distancias intergénicas.
    
    Args:
        distancias: Lista de distancias intergénicas
    
    Returns:
        dict: Conteo de cada tipo
    """
    tipos = [d["tipo"] for d in distancias]
    distribucion = Counter(tipos)
    
    print("\n" + "=" * 60)
    print("DISTRIBUCIÓN DE TIPOS DE DISTANCIAS")
    print("=" * 60)
    
    total = len(distancias)
    
    for tipo, conteo in distribucion.most_common():
        porcentaje = (conteo / total) * 100
        print(f"{tipo:<25} {conteo:>6,} ({porcentaje:>5.2f}%)")
    
    return dict(distribucion)


# FUNCIONES DE EXPORTACIÓN
# =============================================================================

def exportar_distancias_csv(distancias, archivo_salida):
    """
    Exporta las distancias intergénicas a un archivo CSV.
    
    Args:
        distancias: Lista de distancias intergénicas
        archivo_salida: Ruta del archivo CSV de salida
    """
    print(f"\nExportando a CSV: {os.path.basename(archivo_salida)}")
    
    with open(archivo_salida, 'w', newline='', encoding='utf-8') as f:
        campos = [
            "gen1_locus", "gen1_nombre", "gen1_fin", "gen1_hebra",
            "gen2_locus", "gen2_nombre", "gen2_inicio", "gen2_hebra",
            "distancia_pb", "tipo", "misma_hebra"
        ]
        
        writer = csv.DictWriter(f, fieldnames=campos)
        writer.writeheader()
        writer.writerows(distancias)
    
    print(f"[OK] {len(distancias):,} distancias exportadas a CSV")


def exportar_estadisticas_json(estadisticas_general, estadisticas_hebra, 
                               distribucion_tipos, regiones_grandes,
                               archivo_salida):
    """
    Exporta todas las estadísticas a un archivo JSON.
    
    Args:
        estadisticas_general: Estadísticas generales
        estadisticas_hebra: Estadísticas por hebra
        distribucion_tipos: Distribución de tipos
        regiones_grandes: Lista de regiones grandes
        archivo_salida: Ruta del archivo JSON de salida
    """
    print(f"\nExportando estadísticas a JSON: {os.path.basename(archivo_salida)}")
    
    # Top 20 regiones grandes (resumido)
    top_regiones = []
    for region in regiones_grandes[:20]:
        top_regiones.append({
            "gen1": region["gen1_nombre"] or region["gen1_locus"],
            "gen2": region["gen2_nombre"] or region["gen2_locus"],
            "distancia_pb": region["distancia_pb"],
            "hebra_gen1": region["gen1_hebra"],
            "hebra_gen2": region["gen2_hebra"]
        })
    
    datos_completos = {
        "metadata": {
            "fecha_analisis": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "script": "analisis_distancias_intergenicas.py",
            "organismo": "Escherichia coli K-12 MG1655"
        },
        "estadisticas_generales": estadisticas_general,
        "estadisticas_por_hebra": estadisticas_hebra,
        "distribucion_tipos": distribucion_tipos,
        "top_20_regiones_grandes": top_regiones
    }
    
    with open(archivo_salida, 'w', encoding='utf-8') as f:
        json.dump(datos_completos, f, indent=2, ensure_ascii=False)
    
    print(f"[OK] Estadísticas exportadas a JSON")


def exportar_regiones_grandes_csv(regiones_grandes, archivo_salida):
    """
    Exporta las regiones intergénicas grandes a un archivo CSV separado.
    
    Args:
        regiones_grandes: Lista de regiones grandes
        archivo_salida: Ruta del archivo CSV de salida
    """
    print(f"\nExportando regiones grandes a CSV: {os.path.basename(archivo_salida)}")
    
    with open(archivo_salida, 'w', newline='', encoding='utf-8') as f:
        campos = [
            "ranking", "distancia_pb", 
            "gen1_locus", "gen1_nombre", "gen1_hebra",
            "gen2_locus", "gen2_nombre", "gen2_hebra",
            "misma_hebra"
        ]
        
        writer = csv.DictWriter(f, fieldnames=campos)
        writer.writeheader()
        
        for i, region in enumerate(regiones_grandes, 1):
            fila = {
                "ranking": i,
                "distancia_pb": region["distancia_pb"],
                "gen1_locus": region["gen1_locus"],
                "gen1_nombre": region["gen1_nombre"],
                "gen1_hebra": region["gen1_hebra"],
                "gen2_locus": region["gen2_locus"],
                "gen2_nombre": region["gen2_nombre"],
                "gen2_hebra": region["gen2_hebra"],
                "misma_hebra": region["misma_hebra"]
            }
            writer.writerow(fila)
    
    print(f"[OK] {len(regiones_grandes):,} regiones grandes exportadas a CSV")


# FUNCIÓN PRINCIPAL
# =============================================================================

def main():
    """
    Función principal que ejecuta todo el análisis.
    """
    print("\n" + "=" * 60)
    print("ANÁLISIS DE DISTANCIAS INTERGÉNICAS - E. coli K-12 MG1655")
    print("=" * 60)
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # 1. Extraer genes ordenados
    genes, longitud_genoma = extraer_genes_ordenados(ARCHIVO_GENBANK)
    
    # 2. Calcular distancias intergénicas
    distancias = calcular_distancias_intergenicas(genes)
    
    # 3. Analizar estadísticas generales
    estadisticas_general = analizar_estadisticas_distancias(distancias)
    
    # 4. Analizar por hebra
    estadisticas_hebra = analizar_por_hebra(distancias)
    
    # 5. Analizar distribución de tipos
    distribucion_tipos = analizar_distribucion_tipos(distancias)
    
    # 6. Encontrar regiones grandes
    regiones_grandes = encontrar_regiones_grandes(distancias)
    
    # 7. Exportar resultados
    print("\n" + "=" * 60)
    print("EXPORTANDO RESULTADOS")
    print("=" * 60)
    
    # Archivo CSV con todas las distancias
    archivo_distancias_csv = os.path.join(RUTA_RESULTADOS, 
                                          "distancias_intergenicas.csv")
    exportar_distancias_csv(distancias, archivo_distancias_csv)
    
    # Archivo CSV solo con regiones grandes
    archivo_grandes_csv = os.path.join(RUTA_RESULTADOS, 
                                       "regiones_intergenicas_grandes.csv")
    exportar_regiones_grandes_csv(regiones_grandes, archivo_grandes_csv)
    
    # Archivo JSON con estadísticas completas
    archivo_json = os.path.join(RUTA_RESULTADOS, 
                                "analisis_distancias_completo.json")
    exportar_estadisticas_json(estadisticas_general, estadisticas_hebra,
                              distribucion_tipos, regiones_grandes,
                              archivo_json)
    
    # Resumen final
    print("\n" + "=" * 60)
    print("ANÁLISIS COMPLETADO")
    print("=" * 60)
    print(f"\nArchivos generados en: {RUTA_RESULTADOS}")
    print(f"  • distancias_intergenicas.csv")
    print(f"  • regiones_intergenicas_grandes.csv")
    print(f"  • analisis_distancias_completo.json")
    print("\n" + "=" * 60 + "\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Búsqueda de genomas en NCBI en tiempo real

Este script busca genomas en la base de datos Nucleotide de NCBI
usando la API Entrez de BioPython y retorna resultados en formato JSON
para consumo de la API web.

Uso:
    python buscar_ncbi.py "Escherichia coli" 20
    python buscar_ncbi.py "Salmonella" 10

Fecha: 2026
Proyecto: GenomeHub - Aplicación web de análisis genómico
"""

from Bio import Entrez
import sys
import json
import urllib.error

# CONFIGURACIÓN
# =============================================================================

# Email requerido por NCBI (obligatorio para usar Entrez)
CORREO_ELECTRONICO = "194522@unsaac.edu.pe"
Entrez.email = CORREO_ELECTRONICO

# FUNCIONES
# =============================================================================

def buscar_genomas(query, limit=20):
    """
    Busca genomas completos en NCBI Nucleotide database.

    Args:
        query (str): Término de búsqueda (nombre organismo, accession, etc.)
        limit (int): Número máximo de resultados (default: 20)

    Returns:
        dict: JSON con resultados estructurados:
            {
                "success": true/false,
                "query": "término buscado",
                "count": N,
                "results": [...]
            }
    """
    try:
        # Búsqueda en Nucleotide database
        # Filtro: solo genomas completos (complete genome)
        search_term = f"{query}[Organism] AND complete genome[Title]"

        search_handle = Entrez.esearch(
            db="nucleotide",
            term=search_term,
            retmax=limit,
            sort="relevance"  # Ordenar por relevancia
        )

        search_results = Entrez.read(search_handle)
        search_handle.close()

        # Lista de IDs encontrados
        ids = search_results["IdList"]

        if not ids:
            return {
                "success": True,
                "query": query,
                "count": 0,
                "results": [],
                "message": "No se encontraron genomas para esta búsqueda"
            }

        # Obtener detalles de los genomas encontrados (summaries)
        # Esto retorna metadata sin descargar los genomas completos
        fetch_handle = Entrez.efetch(
            db="nucleotide",
            id=ids,
            rettype="docsum",  # Document summary (metadata ligera)
            retmode="xml"
        )

        summaries = Entrez.read(fetch_handle)
        fetch_handle.close()

        # Procesar cada resultado
        results = []
        for summary in summaries:
            # Extraer información relevante
            title = summary.get("Title", "")
            organism = extract_organism(title)
            length = int(summary.get("Length", 0))

            result = {
                "id": summary.get("Id", ""),
                "accession": summary.get("AccessionVersion", "N/A"),
                "title": title,
                "organism": organism,
                "length": length,
                "length_mb": round(length / 1e6, 2),
                "update_date": summary.get("UpdateDate", "N/A"),
                "create_date": summary.get("CreateDate", "N/A")
            }

            results.append(result)

        return {
            "success": True,
            "query": query,
            "count": len(results),
            "results": results
        }

    except urllib.error.HTTPError as e:
        return {
            "success": False,
            "error": f"Error HTTP al conectar con NCBI: {e.code}",
            "query": query
        }

    except urllib.error.URLError as e:
        return {
            "success": False,
            "error": f"Error de conexión con NCBI: {e.reason}",
            "query": query
        }

    except Exception as e:
        return {
            "success": False,
            "error": f"Error inesperado: {str(e)}",
            "query": query
        }


def extract_organism(title):
    """
    Extrae el nombre del organismo del título del genoma.

    Los títulos de NCBI suelen tener el formato:
    "Organismo nombre strain, complete genome"

    Args:
        title (str): Título del genoma

    Returns:
        str: Nombre del organismo extraído
    """
    # Típicamente el organismo está antes de la primera coma
    if "," in title:
        organism = title.split(",")[0].strip()
    else:
        # Si no hay coma, tomar primeras 50 caracteres
        organism = title[:50].strip()

    # Remover frases comunes al final
    organism = organism.replace(" complete genome", "")
    organism = organism.replace(" complete sequence", "")

    return organism


# MAIN
# =============================================================================

if __name__ == "__main__":
    # Validar argumentos
    if len(sys.argv) < 2:
        error_result = {
            "success": False,
            "error": "Falta el parámetro de búsqueda (query)",
            "usage": "python buscar_ncbi.py <query> [limit]"
        }
        print(json.dumps(error_result, ensure_ascii=False, indent=2))
        sys.exit(1)

    # Parámetros
    query = sys.argv[1]
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 20

    # Validar limit
    if limit < 1:
        limit = 1
    if limit > 100:
        limit = 100

    # Ejecutar búsqueda
    result = buscar_genomas(query, limit)

    # Retornar JSON (stdout para que PHP lo capture)
    print(json.dumps(result, ensure_ascii=False, indent=2))

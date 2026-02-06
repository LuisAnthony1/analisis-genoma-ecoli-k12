#!/usr/bin/env python3
"""
Descarga de genomas desde NCBI (versión API no-interactiva)

Versión modificada de descargar_genoma.py para uso desde la API web.
Recibe parámetros por línea de comandos y retorna JSON.

Uso:
    python descargar_genoma_api.py "NC_000913.3" "Escherichia coli K-12"
    python descargar_genoma_api.py "NC_003197.2" "Salmonella enterica LT2"

Fecha: 2026
Proyecto: GenomeHub - Aplicación web de análisis genómico
"""

from Bio import Entrez, SeqIO
import os
import sys
import json
import re
from datetime import datetime
import urllib.error

# CONFIGURACIÓN
# =============================================================================

# Email requerido por NCBI
CORREO_ELECTRONICO = "194522@unsaac.edu.pe"
Entrez.email = CORREO_ELECTRONICO

# Rutas de archivos
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_BASE, "crudo")

# FUNCIONES
# =============================================================================

def crear_directorios():
    """Crea el directorio de datos si no existe."""
    if not os.path.exists(RUTA_DATOS_CRUDO):
        os.makedirs(RUTA_DATOS_CRUDO)


def sanitize_filename(name):
    """
    Crea un nombre de archivo seguro desde el nombre del organismo.

    Args:
        name (str): Nombre del organismo

    Returns:
        str: Nombre de archivo sanitizado
    """
    # Convertir a minúsculas
    safe_name = name.lower()

    # Reemplazar espacios y caracteres especiales por guión bajo
    safe_name = re.sub(r'[^a-z0-9]+', '_', safe_name)

    # Remover guiones bajos al inicio/final
    safe_name = safe_name.strip('_')

    # Limitar longitud
    if len(safe_name) > 50:
        safe_name = safe_name[:50]

    return safe_name


def descargar_genoma(accession_id, organism_name):
    """
    Descarga un genoma desde NCBI por accession ID.

    Args:
        accession_id (str): ID de acceso NCBI (ej: "NC_000913.3")
        organism_name (str): Nombre del organismo

    Returns:
        dict: JSON con resultado de la descarga
    """
    try:
        # Crear nombre de archivo seguro
        safe_name = sanitize_filename(organism_name)
        archivo_gb = os.path.join(RUTA_DATOS_CRUDO, f"{safe_name}.gb")
        archivo_fasta = os.path.join(RUTA_DATOS_CRUDO, f"{safe_name}.fasta")
        archivo_metadata = os.path.join(RUTA_DATOS_CRUDO, f"metadata_{safe_name}.json")

        # Crear directorios si no existen
        crear_directorios()

        # =========================================================================
        # PASO 1: Descargar GenBank (con anotaciones completas)
        # =========================================================================

        # Usar "gbwithparts" para genomas grandes con muchas anotaciones
        handle_gb = Entrez.efetch(
            db="nucleotide",
            id=accession_id,
            rettype="gbwithparts",
            retmode="text"
        )

        contenido_gb = handle_gb.read()
        handle_gb.close()

        # Verificar que tiene anotaciones CDS
        if "CDS" not in contenido_gb:
            # Intentar formato "gb" alternativo
            handle_gb_alt = Entrez.efetch(
                db="nucleotide",
                id=accession_id,
                rettype="gb",
                retmode="text"
            )
            contenido_gb_alt = handle_gb_alt.read()
            handle_gb_alt.close()

            if "CDS" in contenido_gb_alt:
                contenido_gb = contenido_gb_alt

        # Guardar GenBank
        with open(archivo_gb, "w") as f:
            f.write(contenido_gb)

        # =========================================================================
        # PASO 2: Descargar FASTA (secuencia pura)
        # =========================================================================

        handle_fasta = Entrez.efetch(
            db="nucleotide",
            id=accession_id,
            rettype="fasta",
            retmode="text"
        )

        contenido_fasta = handle_fasta.read()
        handle_fasta.close()

        # Guardar FASTA
        with open(archivo_fasta, "w") as f:
            f.write(contenido_fasta)

        # =========================================================================
        # PASO 3: Validar y extraer metadata
        # =========================================================================

        # Leer el archivo GenBank para extraer información
        registro = SeqIO.read(archivo_gb, "genbank")

        informacion_genoma = {
            "id": registro.id,
            "nombre": registro.name,
            "descripcion": registro.description,
            "longitud": len(registro.seq),
            "organismo": registro.annotations.get("organism", organism_name),
            "fecha_actualizacion": registro.annotations.get("date", ""),
            "num_features": len(registro.features),
            "nombre_corto": safe_name
        }

        # =========================================================================
        # PASO 4: Guardar metadata JSON
        # =========================================================================

        metadata = {
            "fecha_descarga": datetime.now().isoformat(),
            "fuente": "NCBI",
            "id_acceso": accession_id,
            "email_usado": CORREO_ELECTRONICO,
            "organismo_config": organism_name,
            "genoma": informacion_genoma
        }

        with open(archivo_metadata, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)

        # =========================================================================
        # RETORNAR RESULTADO EXITOSO
        # =========================================================================

        tamano_gb_mb = os.path.getsize(archivo_gb) / (1024 * 1024)
        tamano_fasta_mb = os.path.getsize(archivo_fasta) / (1024 * 1024)

        return {
            "success": True,
            "message": f"Genoma descargado exitosamente: {organism_name}",
            "accession_id": accession_id,
            "organism": organism_name,
            "files": {
                "genbank": archivo_gb,
                "fasta": archivo_fasta,
                "metadata": archivo_metadata
            },
            "sizes": {
                "genbank_mb": round(tamano_gb_mb, 2),
                "fasta_mb": round(tamano_fasta_mb, 2)
            },
            "metadata": metadata
        }

    except urllib.error.HTTPError as e:
        return {
            "success": False,
            "error": f"Error HTTP {e.code}: No se pudo descargar desde NCBI",
            "details": f"Accession ID: {accession_id} - Verifica que sea válido"
        }

    except urllib.error.URLError as e:
        return {
            "success": False,
            "error": f"Error de conexión: {e.reason}",
            "details": "Verifica tu conexión a internet"
        }

    except FileNotFoundError as e:
        return {
            "success": False,
            "error": "Error al parsear el archivo GenBank descargado",
            "details": str(e)
        }

    except Exception as e:
        return {
            "success": False,
            "error": f"Error inesperado: {str(e)}",
            "details": f"Tipo de error: {type(e).__name__}"
        }


# MAIN
# =============================================================================

if __name__ == "__main__":
    # Validar argumentos
    if len(sys.argv) < 3:
        error_result = {
            "success": False,
            "error": "Faltan parámetros requeridos",
            "usage": "python descargar_genoma_api.py <accession_id> <organism_name>",
            "example": "python descargar_genoma_api.py NC_000913.3 'Escherichia coli K-12'"
        }
        print(json.dumps(error_result, ensure_ascii=False, indent=2))
        sys.exit(1)

    # Parámetros
    accession_id = sys.argv[1]
    organism_name = sys.argv[2]

    # Ejecutar descarga
    result = descargar_genoma(accession_id, organism_name)

    # Retornar JSON (stdout para que PHP lo capture)
    print(json.dumps(result, ensure_ascii=False, indent=2))

    # Exit code según resultado
    sys.exit(0 if result["success"] else 1)

#!/usr/bin/env python3
"""
Descarga del Genoma de E. coli K-12 MG1655 desde NCBI

Este script descarga programaticamente el genoma completo de E. coli K-12 MG1655
utilizando las APIs de NCBI a traves del modulo Entrez de BioPython.

Fecha: 2026
Proyecto: Analisis del genoma de E. coli K-12 MG1655
"""

from Bio import Entrez, SeqIO
import os
import sys
import json
from datetime import datetime


# CONFIGURACION
# =============================================================================

# Configuracion de NCBI Entrez (IMPORTANTE: usar tu email real)
CORREO_ELECTRONICO = "194522@unsaac.edu.pe"
Entrez.email = CORREO_ELECTRONICO

# ID del genoma de E. coli K-12 MG1655 en NCBI
ID_GENOMA = "U00096.3"  # Accession number en GenBank

# Rutas de archivos
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_BASE, "datos", "crudo")
ARCHIVO_GENBANK = os.path.join(RUTA_DATOS_CRUDO, "ecoli_k12.gb")
ARCHIVO_FASTA = os.path.join(RUTA_DATOS_CRUDO, "ecoli_k12.fasta")
ARCHIVO_METADATA = os.path.join(RUTA_DATOS_CRUDO, "metadata_descarga.json")

# Valores esperados para validacion
VALORES_ESPERADOS = {
    "organismo": "Escherichia coli",
    "cepa": "K-12",
    "longitud_aproximada": 4641652,  # pares de bases
    "tolerancia_longitud": 1000  # margen de error permitido
}


# FUNCIONES
# =============================================================================

def crear_directorios():
    """
    Crea los directorios necesarios si no existen.
    """
    if not os.path.exists(RUTA_DATOS_CRUDO):
        os.makedirs(RUTA_DATOS_CRUDO)
        print(f"[INFO] Directorio creado: {RUTA_DATOS_CRUDO}")


def descargar_genoma(formato="gb"):
    """
    Descarga el genoma desde NCBI en el formato especificado.

    Args:
        formato: "gb" para GenBank o "fasta" para FASTA

    Returns:
        str: Contenido del archivo descargado

    Raises:
        Exception: Si hay errores de conectividad o el genoma no se encuentra
    """
    nombre_formato = "GenBank (incluye anotaciones)" if formato == "gb" else "FASTA (secuencia pura)"
    print(f"[INFO] Descargando genoma {ID_GENOMA} en formato {nombre_formato}...")

    try:
        # Realizar la peticion a NCBI
        manejador = Entrez.efetch(
            db="nucleotide",
            id=ID_GENOMA,
            rettype=formato,
            retmode="text"
        )

        contenido = manejador.read()
        manejador.close()

        print(f"[OK] Descarga {formato.upper()} completada")
        return contenido

    except Exception as error:
        print(f"[ERROR] Error al descargar el genoma ({formato}): {error}")
        raise


def guardar_archivo(contenido, ruta_archivo):
    """
    Guarda el contenido descargado en un archivo.

    Args:
        contenido: Contenido a guardar
        ruta_archivo: Ruta donde guardar el archivo
    """
    print(f"[INFO] Guardando archivo en: {ruta_archivo}")

    with open(ruta_archivo, "w") as archivo:
        archivo.write(contenido)

    tamano_mb = os.path.getsize(ruta_archivo) / (1024 * 1024)
    print(f"[OK] Archivo guardado ({tamano_mb:.2f} MB)")


def validar_genoma(ruta_archivo):
    """
    Valida que el genoma descargado corresponde a E. coli K-12 MG1655.

    Verifica:
    - Que el organismo sea Escherichia coli
    - Que la cepa sea K-12
    - Que la longitud sea aproximadamente 4.6 millones de pares de bases

    Args:
        ruta_archivo: Ruta al archivo GenBank descargado

    Returns:
        dict: Informacion del genoma validado
    """
    print("[INFO] Validando integridad del genoma descargado...")

    # Leer el archivo GenBank
    registro = SeqIO.read(ruta_archivo, "genbank")

    # Extraer informacion
    informacion = {
        "id": registro.id,
        "nombre": registro.name,
        "descripcion": registro.description,
        "longitud": len(registro.seq),
        "organismo": registro.annotations.get("organism", "Desconocido"),
        "fecha_actualizacion": registro.annotations.get("date", "Desconocida"),
        "num_features": len(registro.features)
    }

    # Validar organismo
    if VALORES_ESPERADOS["organismo"] not in informacion["organismo"]:
        print(f"[ADVERTENCIA] Organismo inesperado: {informacion['organismo']}")
    else:
        print(f"[OK] Organismo correcto: {informacion['organismo']}")

    # Validar longitud
    diferencia = abs(informacion["longitud"] - VALORES_ESPERADOS["longitud_aproximada"])
    if diferencia > VALORES_ESPERADOS["tolerancia_longitud"]:
        print(f"[ADVERTENCIA] Longitud inesperada: {informacion['longitud']:,} pb")
        print(f"    Esperado: ~{VALORES_ESPERADOS['longitud_aproximada']:,} pb")
    else:
        print(f"[OK] Longitud correcta: {informacion['longitud']:,} pares de bases")

    # Validar que tenga anotaciones
    if informacion["num_features"] > 0:
        print(f"[OK] Anotaciones encontradas: {informacion['num_features']:,} features")
    else:
        print("[ADVERTENCIA] No se encontraron anotaciones en el archivo")

    return informacion


def guardar_metadata(informacion):
    """
    Guarda los metadatos de la descarga en formato JSON.

    Args:
        informacion: Diccionario con informacion del genoma
    """
    metadata = {
        "fecha_descarga": datetime.now().isoformat(),
        "fuente": "NCBI",
        "id_acceso": ID_GENOMA,
        "email_usado": CORREO_ELECTRONICO,
        "genoma": informacion
    }

    with open(ARCHIVO_METADATA, "w", encoding="utf-8") as archivo:
        json.dump(metadata, archivo, indent=2, ensure_ascii=False)

    print(f"[OK] Metadata guardada en: {ARCHIVO_METADATA}")


def mostrar_resumen(informacion):
    """
    Muestra un resumen de la descarga.

    Args:
        informacion: Diccionario con informacion del genoma
    """
    print("\n" + "=" * 60)
    print("RESUMEN DE LA DESCARGA")
    print("=" * 60)
    print(f"  ID:          {informacion['id']}")
    print(f"  Organismo:   {informacion['organismo']}")
    print(f"  Descripcion: {informacion['descripcion'][:50]}...")
    print(f"  Longitud:    {informacion['longitud']:,} pb")
    print(f"  Features:    {informacion['num_features']:,}")
    print(f"  GenBank:     {ARCHIVO_GENBANK}")
    print(f"  FASTA:       {ARCHIVO_FASTA}")
    print("=" * 60)


# EJECUCION PRINCIPAL
# =============================================================================

def main():
    """
    Funcion principal que ejecuta todo el proceso de descarga.
    """
    print("\n" + "=" * 60)
    print("DESCARGA DEL GENOMA DE E. coli K-12 MG1655")
    print("=" * 60 + "\n")

    # Verificar que se haya configurado el email
    if CORREO_ELECTRONICO == "tu_correo@ejemplo.com":
        print("[ERROR] Debes configurar tu correo electronico en CORREO_ELECTRONICO")
        print("        NCBI requiere un email valido para usar sus APIs")
        sys.exit(1)

    # Paso 1: Crear directorios
    crear_directorios()

    # Paso 2: Descargar genoma en formato GenBank
    try:
        contenido_genbank = descargar_genoma(formato="gb")
    except Exception as error:
        print(f"[ERROR] No se pudo completar la descarga GenBank: {error}")
        sys.exit(1)

    # Paso 3: Descargar genoma en formato FASTA
    try:
        contenido_fasta = descargar_genoma(formato="fasta")
    except Exception as error:
        print(f"[ERROR] No se pudo completar la descarga FASTA: {error}")
        sys.exit(1)

    # Paso 4: Guardar archivos
    guardar_archivo(contenido_genbank, ARCHIVO_GENBANK)
    guardar_archivo(contenido_fasta, ARCHIVO_FASTA)

    # Paso 5: Validar genoma
    informacion_genoma = validar_genoma(ARCHIVO_GENBANK)

    # Paso 6: Guardar metadata
    guardar_metadata(informacion_genoma)

    # Paso 7: Mostrar resumen
    mostrar_resumen(informacion_genoma)

    print("\n[OK] Proceso completado exitosamente\n")


if __name__ == "__main__":
    main()

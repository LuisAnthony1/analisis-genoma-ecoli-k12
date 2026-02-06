#!/usr/bin/env python3
"""
Descarga de Genomas desde NCBI

Este script descarga programaticamente genomas completos utilizando
las APIs de NCBI a traves del modulo Entrez de BioPython.

Organismos disponibles:
1. Escherichia coli K-12 MG1655 (cepa de laboratorio de referencia)
2. Salmonella enterica serovar Typhimurium LT2 (para comparacion)

Fecha: 2026
Proyecto: Analisis comparativo de genomas bacterianos
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

# Rutas de archivos
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_BASE, "datos", "crudo")

# Configuracion de organismos disponibles
ORGANISMOS = {
    1: {
        "nombre": "Escherichia coli K-12 MG1655",
        "nombre_corto": "ecoli_k12",
        "id_genoma": "U00096.3",  # Accession number en GenBank
        "descripcion": "Cepa de laboratorio de referencia, no patogena",
        "longitud_esperada": 4641652,
        "tolerancia": 1000,
        "archivo_gb": "ecoli_k12.gb",
        "archivo_fasta": "ecoli_k12.fasta",
        "archivo_metadata": "metadata_ecoli.json"
    },
    2: {
        "nombre": "Salmonella enterica serovar Typhimurium LT2",
        "nombre_corto": "salmonella_lt2",
        "id_genoma": "NC_003197.2",  # Genoma de referencia NCBI
        "descripcion": "Patogeno - causa gastroenteritis, cepa de referencia",
        "longitud_esperada": 4857450,  # ~4.86 Mb
        "tolerancia": 5000,
        "archivo_gb": "salmonella_lt2.gb",
        "archivo_fasta": "salmonella_lt2.fasta",
        "archivo_metadata": "metadata_salmonella.json"
    }
}


# FUNCIONES
# =============================================================================

def mostrar_menu():
    """
    Muestra el menu de seleccion de organismos.

    Returns:
        int: Opcion seleccionada (1, 2, o 3 para ambos)
    """
    print("\n" + "=" * 70)
    print("DESCARGA DE GENOMAS BACTERIANOS DESDE NCBI")
    print("=" * 70)
    print("\n  Seleccione el genoma a descargar:\n")

    for num, org in ORGANISMOS.items():
        print(f"    [{num}] {org['nombre']}")
        print(f"        {org['descripcion']}")
        print(f"        ID: {org['id_genoma']} | Tamano: ~{org['longitud_esperada']/1e6:.2f} Mb")
        print()

    print(f"    [3] Descargar AMBOS genomas (para comparacion)")
    print(f"    [0] Salir")
    print()

    while True:
        try:
            opcion = int(input("  Ingrese su opcion (0-3): "))
            if opcion in [0, 1, 2, 3]:
                return opcion
            print("  [ERROR] Opcion no valida. Ingrese 0, 1, 2 o 3.")
        except ValueError:
            print("  [ERROR] Ingrese un numero valido.")


def crear_directorios():
    """
    Crea los directorios necesarios si no existen.
    """
    if not os.path.exists(RUTA_DATOS_CRUDO):
        os.makedirs(RUTA_DATOS_CRUDO)
        print(f"[INFO] Directorio creado: {RUTA_DATOS_CRUDO}")


def descargar_genoma(id_genoma, nombre_organismo, formato="gb"):
    """
    Descarga el genoma desde NCBI en el formato especificado.

    Args:
        id_genoma: ID de acceso del genoma en NCBI
        nombre_organismo: Nombre del organismo (para mensajes)
        formato: "gb" para GenBank o "fasta" para FASTA

    Returns:
        str: Contenido del archivo descargado

    Raises:
        Exception: Si hay errores de conectividad o el genoma no se encuentra
    """
    nombre_formato = "GenBank (incluye anotaciones)" if formato == "gb" else "FASTA (secuencia pura)"
    print(f"\n[INFO] Descargando {nombre_organismo}")
    print(f"       ID: {id_genoma} | Formato: {nombre_formato}")

    try:
        # Para GenBank, usar "gbwithparts" para obtener todas las anotaciones
        # Esto es necesario para genomas grandes con muchas features
        if formato == "gb":
            rettype = "gbwithparts"
        else:
            rettype = formato

        # Realizar la peticion a NCBI
        manejador = Entrez.efetch(
            db="nucleotide",
            id=id_genoma,
            rettype=rettype,
            retmode="text"
        )

        contenido = manejador.read()
        manejador.close()

        # Verificar que el GenBank tiene anotaciones
        if formato == "gb" and "CDS" not in contenido:
            print(f"[WARN] El archivo descargado puede no tener anotaciones completas")
            print(f"       Intentando descarga alternativa...")

            # Intentar con formato "gb" directamente
            manejador = Entrez.efetch(
                db="nucleotide",
                id=id_genoma,
                rettype="gb",
                retmode="text"
            )
            contenido_alt = manejador.read()
            manejador.close()

            if "CDS" in contenido_alt:
                contenido = contenido_alt
                print(f"[OK] Descarga alternativa exitosa con anotaciones")

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


def validar_genoma(ruta_archivo, config_organismo):
    """
    Valida que el genoma descargado corresponde al organismo esperado.

    Args:
        ruta_archivo: Ruta al archivo GenBank descargado
        config_organismo: Configuracion del organismo

    Returns:
        dict: Informacion del genoma validado
    """
    print(f"\n[INFO] Validando integridad del genoma...")

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
        "num_features": len(registro.features),
        "nombre_corto": config_organismo["nombre_corto"]
    }

    # Validar longitud
    diferencia = abs(informacion["longitud"] - config_organismo["longitud_esperada"])
    if diferencia > config_organismo["tolerancia"]:
        print(f"[ADVERTENCIA] Longitud inesperada: {informacion['longitud']:,} pb")
        print(f"    Esperado: ~{config_organismo['longitud_esperada']:,} pb")
    else:
        print(f"[OK] Organismo: {informacion['organismo']}")
        print(f"[OK] Longitud: {informacion['longitud']:,} pares de bases")
        print(f"[OK] Features: {informacion['num_features']:,} anotaciones")

    return informacion


def guardar_metadata(informacion, ruta_archivo, config_organismo):
    """
    Guarda los metadatos de la descarga en formato JSON.

    Args:
        informacion: Diccionario con informacion del genoma
        ruta_archivo: Ruta donde guardar el JSON
        config_organismo: Configuracion del organismo
    """
    metadata = {
        "fecha_descarga": datetime.now().isoformat(),
        "fuente": "NCBI",
        "id_acceso": config_organismo["id_genoma"],
        "email_usado": CORREO_ELECTRONICO,
        "organismo_config": config_organismo["nombre"],
        "genoma": informacion
    }

    with open(ruta_archivo, "w", encoding="utf-8") as archivo:
        json.dump(metadata, archivo, indent=2, ensure_ascii=False)

    print(f"[OK] Metadata guardada: {os.path.basename(ruta_archivo)}")


def descargar_organismo(num_organismo):
    """
    Descarga un organismo completo (GenBank + FASTA + validacion).

    Args:
        num_organismo: Numero del organismo (1 o 2)

    Returns:
        dict: Informacion del genoma descargado
    """
    config = ORGANISMOS[num_organismo]

    print("\n" + "=" * 70)
    print(f"DESCARGANDO: {config['nombre']}")
    print("=" * 70)

    # Rutas de archivos
    archivo_gb = os.path.join(RUTA_DATOS_CRUDO, config["archivo_gb"])
    archivo_fasta = os.path.join(RUTA_DATOS_CRUDO, config["archivo_fasta"])
    archivo_metadata = os.path.join(RUTA_DATOS_CRUDO, config["archivo_metadata"])

    # Descargar GenBank
    try:
        contenido_genbank = descargar_genoma(
            config["id_genoma"],
            config["nombre"],
            formato="gb"
        )
    except Exception as error:
        print(f"[ERROR] No se pudo descargar GenBank: {error}")
        return None

    # Descargar FASTA
    try:
        contenido_fasta = descargar_genoma(
            config["id_genoma"],
            config["nombre"],
            formato="fasta"
        )
    except Exception as error:
        print(f"[ERROR] No se pudo descargar FASTA: {error}")
        return None

    # Guardar archivos
    guardar_archivo(contenido_genbank, archivo_gb)
    guardar_archivo(contenido_fasta, archivo_fasta)

    # Validar genoma
    informacion = validar_genoma(archivo_gb, config)

    # Guardar metadata
    guardar_metadata(informacion, archivo_metadata, config)

    return informacion


def mostrar_resumen(organismos_descargados):
    """
    Muestra un resumen de las descargas completadas.

    Args:
        organismos_descargados: Lista de diccionarios con informacion
    """
    print("\n" + "=" * 70)
    print("RESUMEN DE DESCARGAS COMPLETADAS")
    print("=" * 70)

    for info in organismos_descargados:
        if info:
            print(f"\n  {info['organismo']}")
            print(f"  " + "-" * 50)
            print(f"    ID:        {info['id']}")
            print(f"    Longitud:  {info['longitud']:,} pb ({info['longitud']/1e6:.2f} Mb)")
            print(f"    Features:  {info['num_features']:,}")

    print("\n" + "=" * 70)

    if len(organismos_descargados) == 2:
        print("\n  [INFO] Ambos genomas descargados. Puedes ejecutar:")
        print("         - python analisis_genes.py    (para analizar cada genoma)")
        print("         - python comparar_genomas.py  (para comparar ambos)")

    print()


# EJECUCION PRINCIPAL
# =============================================================================

def main():
    """
    Funcion principal que ejecuta todo el proceso de descarga.
    """
    # Verificar que se haya configurado el email
    if CORREO_ELECTRONICO == "tu_correo@ejemplo.com":
        print("[ERROR] Debes configurar tu correo electronico en CORREO_ELECTRONICO")
        print("        NCBI requiere un email valido para usar sus APIs")
        sys.exit(1)

    # Mostrar menu y obtener opcion
    opcion = mostrar_menu()

    if opcion == 0:
        print("\n[INFO] Saliendo...\n")
        sys.exit(0)

    # Crear directorios
    crear_directorios()

    # Descargar segun opcion
    organismos_descargados = []

    if opcion == 1:
        # Solo E. coli
        info = descargar_organismo(1)
        if info:
            organismos_descargados.append(info)

    elif opcion == 2:
        # Solo Salmonella
        info = descargar_organismo(2)
        if info:
            organismos_descargados.append(info)

    elif opcion == 3:
        # Ambos
        print("\n[INFO] Descargando ambos genomas para comparacion...")

        info_ecoli = descargar_organismo(1)
        if info_ecoli:
            organismos_descargados.append(info_ecoli)

        info_salmonella = descargar_organismo(2)
        if info_salmonella:
            organismos_descargados.append(info_salmonella)

    # Mostrar resumen
    if organismos_descargados:
        mostrar_resumen(organismos_descargados)
        print("[OK] Proceso completado exitosamente\n")
    else:
        print("\n[ERROR] No se pudieron descargar los genomas\n")
        sys.exit(1)


if __name__ == "__main__":
    main()

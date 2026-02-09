#!/usr/bin/env python3
"""
GenomeHub - Consulta de Literatura via IA (Gemini)

Obtiene valores de referencia genomica para cualquier organismo
usando Gemini API, y los guarda como JSON.

USO:
    python consultar_literatura_ia.py <nombre_organismo> [basename]
"""

import json
import os
import re
import sys

# Rutas
DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # backend/
DIRECTORIO_PROYECTO = os.path.dirname(DIRECTORIO_BASE)  # raiz

# Agregar raiz al path para importar gemini_client
sys.path.insert(0, DIRECTORIO_PROYECTO)
from backend.gemini_client import consultar_gemini

RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados", "tablas")


def obtener_literatura(nombre_organismo, basename=None):
    """Consulta Gemini para obtener valores de referencia genomica."""

    prompt = f"""Eres un experto en genomica bacteriana. Necesito los valores de referencia
de la literatura cientifica para el organismo: {nombre_organismo}

Responde UNICAMENTE con un JSON valido (sin markdown, sin explicaciones) con esta estructura exacta:
{{
    "organismo": "{nombre_organismo}",
    "genes_codificantes": <numero entero de genes codificantes conocidos>,
    "gc_porcentaje": <porcentaje GC del genoma, decimal con 1 decimal>,
    "tamano_genoma_pb": <tamano del genoma en pares de bases>,
    "densidad_genica_porcentaje": <porcentaje del genoma que codifica proteinas, decimal>,
    "tamano_promedio_gen_pb": <tamano promedio de gen en pares de bases>,
    "genes_rrna": <numero de genes rRNA>,
    "genes_trna": <numero de genes tRNA>,
    "porcentaje_solapamiento": <porcentaje de genes solapados>,
    "fuente": "Valores de referencia obtenidos via IA basados en literatura cientifica",
    "referencias": [
        "Referencia bibliografica 1",
        "Referencia bibliografica 2"
    ]
}}

IMPORTANTE: Usa datos reales de la literatura. Si no conoces un valor exacto, usa el valor mas cercano conocido.
Responde SOLO con el JSON, sin texto adicional."""

    print(f"[LITERATURA-IA] Consultando Gemini para: {nombre_organismo}...")
    respuesta, error = consultar_gemini(prompt)

    if error:
        print(f"[LITERATURA-IA] Error: {error}")
        return None

    if not respuesta:
        print("[LITERATURA-IA] Respuesta vacia")
        return None

    # Limpiar respuesta (puede venir con ```json ... ```)
    texto = respuesta.strip()
    if texto.startswith("```"):
        texto = re.sub(r'^```\w*\n?', '', texto)
        texto = re.sub(r'\n?```$', '', texto)
        texto = texto.strip()

    try:
        datos = json.loads(texto)
        print(f"[LITERATURA-IA] Datos obtenidos exitosamente")
        return datos
    except json.JSONDecodeError as e:
        print(f"[LITERATURA-IA] Error parseando JSON: {e}")
        print(f"[LITERATURA-IA] Respuesta raw: {texto[:200]}")
        return None


def main():
    if len(sys.argv) < 2:
        print("USO: python consultar_literatura_ia.py <nombre_organismo> [basename]")
        sys.exit(1)

    nombre_organismo = sys.argv[1]
    basename = sys.argv[2] if len(sys.argv) > 2 else None

    if not basename:
        basename = re.sub(r'[^a-z0-9]+', '_', nombre_organismo.lower()).strip('_')[:50]

    datos = obtener_literatura(nombre_organismo, basename)

    if datos is None:
        print("[LITERATURA-IA] No se pudieron obtener datos de literatura")
        sys.exit(1)

    # Guardar resultado
    os.makedirs(RUTA_RESULTADOS, exist_ok=True)
    archivo_salida = os.path.join(RUTA_RESULTADOS, f"literatura_{basename}.json")

    with open(archivo_salida, "w", encoding="utf-8") as f:
        json.dump(datos, f, indent=2, ensure_ascii=False)

    print(f"[LITERATURA-IA] Guardado en: {archivo_salida}")
    print(json.dumps(datos, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

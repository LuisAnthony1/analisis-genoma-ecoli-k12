"""
GenomeHub - Cliente Gemini API con rotacion de claves

Lee claves desde claves_api.txt, selecciona aleatoriamente,
y hace fallback a otra clave si hay error de cuota (429).
"""

import json
import os
import random
import time
import urllib.request
import urllib.error

DIRECTORIO_PROYECTO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARCHIVO_CLAVES = os.path.join(DIRECTORIO_PROYECTO, "claves_api.txt")
GEMINI_MODEL = "gemini-2.5-flash-lite"
GEMINI_URL = f"https://generativelanguage.googleapis.com/v1beta/models/{GEMINI_MODEL}:generateContent"

# Throttle: minimo 1 segundo entre llamadas para no sobrecargar la API
_THROTTLE_SEGUNDOS = 1.0
_ultimo_llamado = 0.0


def cargar_claves():
    """Lee claves API desde claves_api.txt, ignora comentarios y lineas vacias."""
    if not os.path.exists(ARCHIVO_CLAVES):
        return []

    claves = []
    with open(ARCHIVO_CLAVES, "r", encoding="utf-8") as f:
        for linea in f:
            linea = linea.strip()
            if linea and not linea.startswith("#"):
                claves.append(linea)
    return claves


def consultar_gemini(prompt, contexto=""):
    """
    Envia prompt a Gemini API. Selecciona clave aleatoria.
    Si falla por cuota (429), prueba con otra clave.
    Incluye throttle automatico de 1s entre llamadas.

    Returns:
        str: Texto de respuesta o None si falla
        str: Error message si falla, None si exito
    """
    global _ultimo_llamado

    # Throttle: esperar si la ultima llamada fue hace menos de _THROTTLE_SEGUNDOS
    ahora = time.time()
    transcurrido = ahora - _ultimo_llamado
    if transcurrido < _THROTTLE_SEGUNDOS:
        espera = _THROTTLE_SEGUNDOS - transcurrido
        time.sleep(espera)

    _ultimo_llamado = time.time()

    claves = cargar_claves()
    if not claves:
        return None, "No hay claves API configuradas en claves_api.txt"

    mensaje_completo = prompt
    if contexto:
        mensaje_completo = f"Contexto:\n{contexto}\n\nPregunta:\n{prompt}"

    payload = json.dumps({
        "contents": [{
            "parts": [{"text": mensaje_completo}]
        }],
        "generationConfig": {
            "temperature": 0.7,
            "maxOutputTokens": 4096
        }
    }).encode("utf-8")

    # Mezclar claves para seleccion aleatoria
    claves_shuffled = claves[:]
    random.shuffle(claves_shuffled)

    ultimo_error = ""
    for clave in claves_shuffled:
        url = f"{GEMINI_URL}?key={clave}"
        req = urllib.request.Request(
            url,
            data=payload,
            headers={"Content-Type": "application/json"},
            method="POST"
        )

        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                resultado = json.loads(resp.read().decode("utf-8"))

            # Extraer texto de respuesta
            candidates = resultado.get("candidates", [])
            if candidates:
                parts = candidates[0].get("content", {}).get("parts", [])
                if parts:
                    return parts[0].get("text", ""), None

            return None, "Respuesta vacia de Gemini"

        except urllib.error.HTTPError as e:
            body = ""
            try:
                body = e.read().decode("utf-8", errors="replace")[:500]
            except Exception:
                pass
            if e.code == 429:
                ultimo_error = f"Cuota excedida (429) con clave ...{clave[-6:]}"
                continue
            elif e.code == 400:
                return None, f"Error 400: Solicitud invalida - {body}"
            elif e.code == 404:
                return None, f"Error 404: Modelo '{GEMINI_MODEL}' no encontrado - {body}"
            else:
                ultimo_error = f"Error HTTP {e.code}: {body}"
                continue
        except urllib.error.URLError as e:
            return None, f"Sin conexion: {e.reason}"
        except Exception as e:
            ultimo_error = str(e)
            continue

    return None, f"Todas las claves fallaron. Ultimo error: {ultimo_error}"

#!/usr/bin/env python3
"""
GenomeHub - Servidor Web Local
Sirve la aplicación web y conecta con NCBI API para buscar/descargar genomas.

USO:
    python servidor.py

    Luego abre en tu navegador: http://localhost:8080

Requisitos: BioPython (pip install biopython)
"""

import http.server
import json
import os
import sys
import re
import shutil
import subprocess
import time
import urllib.parse
import urllib.error
from datetime import datetime

# Verificar BioPython
try:
    from Bio import Entrez, SeqIO
except ImportError:
    print("ERROR: BioPython no está instalado.")
    print("Instálalo con: pip install biopython")
    sys.exit(1)

# =============================================================================
# CONFIGURACIÓN
# =============================================================================

PUERTO = 8080
CORREO_ELECTRONICO = "194522@unsaac.edu.pe"
Entrez.email = CORREO_ELECTRONICO

# Rutas
DIRECTORIO_BASE = os.path.dirname(os.path.abspath(__file__))
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_BASE, "datos", "crudo")
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "backend", "resultados")

# Directorios legacy donde los scripts escriben por defecto
RUTA_RESULTADOS_TABLAS = os.path.join(RUTA_RESULTADOS, "tablas")
RUTA_RESULTADOS_FIGURAS = os.path.join(RUTA_RESULTADOS, "figuras")


# =============================================================================
# FUNCIONES NCBI
# =============================================================================

def buscar_ncbi(query, limit=20):
    """Busca genomas completos en NCBI."""
    try:
        search_term = f"{query}[Organism] AND complete genome[Title]"

        handle = Entrez.esearch(
            db="nucleotide",
            term=search_term,
            retmax=limit,
            sort="relevance"
        )
        search_results = Entrez.read(handle)
        handle.close()

        ids = search_results["IdList"]

        if not ids:
            return {
                "success": True,
                "query": query,
                "count": 0,
                "results": [],
                "message": "No se encontraron genomas completos para esta búsqueda"
            }

        # Obtener summaries
        fetch_handle = Entrez.efetch(
            db="nucleotide",
            id=ids,
            rettype="docsum",
            retmode="xml"
        )
        summaries = Entrez.read(fetch_handle)
        fetch_handle.close()

        results = []
        for s in summaries:
            title = s.get("Title", "")
            organism = title.split(",")[0].strip() if "," in title else title[:60]
            organism = organism.replace(" complete genome", "").replace(" complete sequence", "")
            length = int(s.get("Length", 0))

            results.append({
                "id": s.get("Id", ""),
                "accession": s.get("AccessionVersion", "N/A"),
                "title": title,
                "organism": organism,
                "length": length,
                "length_mb": round(length / 1e6, 2),
                "update_date": s.get("UpdateDate", "N/A"),
                "create_date": s.get("CreateDate", "N/A")
            })

        return {"success": True, "query": query, "count": len(results), "results": results}

    except urllib.error.HTTPError as e:
        return {"success": False, "error": f"Error HTTP {e.code} al conectar con NCBI"}
    except urllib.error.URLError as e:
        return {"success": False, "error": f"Sin conexión a internet: {e.reason}"}
    except Exception as e:
        return {"success": False, "error": f"Error: {str(e)}"}


def descargar_genoma(accession_id, organism_name):
    """Descarga genoma desde NCBI a datos/crudo/."""
    try:
        # Crear directorio
        os.makedirs(RUTA_DATOS_CRUDO, exist_ok=True)

        # Nombre seguro para archivo
        safe_name = re.sub(r'[^a-z0-9]+', '_', organism_name.lower()).strip('_')[:50]

        archivo_gb = os.path.join(RUTA_DATOS_CRUDO, f"{safe_name}.gb")
        archivo_fasta = os.path.join(RUTA_DATOS_CRUDO, f"{safe_name}.fasta")
        archivo_metadata = os.path.join(RUTA_DATOS_CRUDO, f"metadata_{safe_name}.json")

        # --- DESCARGAR GENBANK ---
        handle = Entrez.efetch(
            db="nucleotide", id=accession_id,
            rettype="gbwithparts", retmode="text"
        )
        contenido_gb = handle.read()
        handle.close()

        # Verificar anotaciones
        if "CDS" not in contenido_gb:
            handle2 = Entrez.efetch(
                db="nucleotide", id=accession_id,
                rettype="gb", retmode="text"
            )
            alt = handle2.read()
            handle2.close()
            if "CDS" in alt:
                contenido_gb = alt

        with open(archivo_gb, "w") as f:
            f.write(contenido_gb)

        # --- DESCARGAR FASTA ---
        handle_f = Entrez.efetch(
            db="nucleotide", id=accession_id,
            rettype="fasta", retmode="text"
        )
        contenido_fasta = handle_f.read()
        handle_f.close()

        with open(archivo_fasta, "w") as f:
            f.write(contenido_fasta)

        # --- EXTRAER METADATA ---
        registro = SeqIO.read(archivo_gb, "genbank")

        info_genoma = {
            "id": registro.id,
            "nombre": registro.name,
            "descripcion": registro.description,
            "longitud": len(registro.seq),
            "organismo": registro.annotations.get("organism", organism_name),
            "fecha_actualizacion": registro.annotations.get("date", ""),
            "num_features": len(registro.features),
            "nombre_corto": safe_name
        }

        metadata = {
            "fecha_descarga": datetime.now().isoformat(),
            "fuente": "NCBI",
            "id_acceso": accession_id,
            "organismo_config": organism_name,
            "genoma": info_genoma
        }

        with open(archivo_metadata, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)

        tam_gb = round(os.path.getsize(archivo_gb) / (1024*1024), 2)
        tam_fasta = round(os.path.getsize(archivo_fasta) / (1024*1024), 2)

        return {
            "success": True,
            "message": f"Genoma descargado: {organism_name}",
            "accession_id": accession_id,
            "organism": organism_name,
            "files": {
                "genbank": f"datos/crudo/{safe_name}.gb",
                "fasta": f"datos/crudo/{safe_name}.fasta",
                "metadata": f"datos/crudo/metadata_{safe_name}.json"
            },
            "sizes": {"genbank_mb": tam_gb, "fasta_mb": tam_fasta},
            "metadata": metadata
        }

    except urllib.error.HTTPError as e:
        return {"success": False, "error": f"Error HTTP {e.code}: Accession ID inválido o no encontrado"}
    except urllib.error.URLError as e:
        return {"success": False, "error": f"Sin conexión a internet: {e.reason}"}
    except Exception as e:
        return {"success": False, "error": f"Error: {str(e)}"}


def listar_genomas():
    """Lista genomas descargados en datos/crudo/."""
    genomas = []

    if not os.path.exists(RUTA_DATOS_CRUDO):
        return {"success": True, "count": 0, "genomes": []}

    for archivo in os.listdir(RUTA_DATOS_CRUDO):
        if not archivo.endswith(".gb"):
            continue

        basename = archivo[:-3]  # sin .gb
        ruta_gb = os.path.join(RUTA_DATOS_CRUDO, archivo)
        ruta_fasta = os.path.join(RUTA_DATOS_CRUDO, f"{basename}.fasta")
        # Buscar metadata en datos/crudo/ o en backend/crudo/
        ruta_meta = None
        posibles_meta = [
            os.path.join(RUTA_DATOS_CRUDO, f"metadata_{basename}.json"),
            os.path.join(DIRECTORIO_BASE, "backend", "crudo", f"metadata_{basename}.json"),
            os.path.join(DIRECTORIO_BASE, "backend", "crudo", f"metadata_{basename.split('_')[0]}.json"),
        ]
        for ruta in posibles_meta:
            if os.path.exists(ruta):
                ruta_meta = ruta
                break

        info = {
            "filename": archivo,
            "basename": basename,
            "size_mb": round(os.path.getsize(ruta_gb) / (1024*1024), 2),
            "has_fasta": os.path.exists(ruta_fasta),
            "has_metadata": ruta_meta is not None
        }

        # Leer metadata si existe
        if ruta_meta:
            try:
                with open(ruta_meta, "r", encoding="utf-8") as f:
                    meta = json.load(f)
                genoma = meta.get("genoma", {})
                info["organism"] = genoma.get("organismo", basename)
                info["accession_id"] = meta.get("id_acceso", "N/A")
                info["length"] = genoma.get("longitud", 0)
                info["length_mb"] = round(genoma.get("longitud", 0) / 1e6, 2)
                info["num_features"] = genoma.get("num_features", 0)
                info["download_date"] = meta.get("fecha_descarga", "")
                info["description"] = genoma.get("descripcion", "")
            except Exception:
                info["organism"] = basename.replace("_", " ").title()
                info["accession_id"] = "N/A"
        else:
            info["organism"] = basename.replace("_", " ").title()
            info["accession_id"] = "N/A"

        genomas.append(info)

    return {"success": True, "count": len(genomas), "genomes": genomas}


def ejecutar_analisis(script_name, organism=None, genome_basename=None, genome_basename_2=None):
    """Ejecuta un script de análisis Python y mueve resultados a carpeta del genoma."""
    scripts_permitidos = {
        "analisis_genes": {"file": "analisis_genes.py", "timeout": 120},
        "analisis_codones": {"file": "analisis_codones.py", "timeout": 60},
        "analisis_distancias_intergenicas": {"file": "analisis_distancias_intergenicas.py", "timeout": 60},
        "comparar_genomas": {"file": "comparar_genomas.py", "timeout": 120},
        "visualizaciones": {"file": "visualizaciones.py", "timeout": 180},
        "consultar_literatura_ia": {"file": "consultar_literatura_ia.py", "timeout": 60},
        "analisis_estructura_gen": {"file": "analisis_estructura_gen.py", "timeout": 180},
        "generar_informe": {"file": "generar_informe.py", "timeout": 300},
    }

    if script_name not in scripts_permitidos:
        return {"success": False, "error": f"Script no permitido. Válidos: {', '.join(scripts_permitidos.keys())}"}

    config = scripts_permitidos[script_name]
    script_path = os.path.join(DIRECTORIO_BASE, "backend", "scripts", config["file"])

    if not os.path.exists(script_path):
        return {"success": False, "error": f"Script no encontrado: {config['file']}"}

    # Construir comando - todos los scripts reciben genome_basename como argumento
    cmd = [sys.executable, script_path]
    if genome_basename:
        cmd.append(genome_basename)
    if genome_basename_2:
        cmd.append(genome_basename_2)

    # Registrar archivos existentes ANTES del análisis
    archivos_antes_tablas = set()
    archivos_antes_figuras = set()
    if os.path.exists(RUTA_RESULTADOS_TABLAS):
        archivos_antes_tablas = set(os.listdir(RUTA_RESULTADOS_TABLAS))
    if os.path.exists(RUTA_RESULTADOS_FIGURAS):
        archivos_antes_figuras = set(os.listdir(RUTA_RESULTADOS_FIGURAS))

    try:
        start = time.time()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=config["timeout"],
            cwd=DIRECTORIO_BASE
        )
        elapsed = round(time.time() - start, 2)

        # Mover archivos nuevos a carpeta del genoma
        if genome_basename:
            mover_resultados_a_genoma(
                genome_basename, archivos_antes_tablas, archivos_antes_figuras
            )

        return {
            "success": True,
            "output": result.stdout + result.stderr,
            "return_code": result.returncode,
            "execution_time": elapsed,
            "script": script_name
        }
    except subprocess.TimeoutExpired:
        return {"success": False, "error": f"Timeout: el script tardó más de {config['timeout']}s"}
    except Exception as e:
        return {"success": False, "error": str(e)}


def mover_resultados_a_genoma(genome_basename, antes_tablas, antes_figuras):
    """Mueve archivos nuevos de resultados a la carpeta del genoma."""
    # Crear directorios destino
    dest_tablas = os.path.join(RUTA_RESULTADOS, genome_basename, "tablas")
    dest_figuras = os.path.join(RUTA_RESULTADOS, genome_basename, "figuras")
    os.makedirs(dest_tablas, exist_ok=True)
    os.makedirs(dest_figuras, exist_ok=True)

    movidos = 0

    # Mover tablas nuevas
    if os.path.exists(RUTA_RESULTADOS_TABLAS):
        for archivo in os.listdir(RUTA_RESULTADOS_TABLAS):
            if archivo not in antes_tablas:
                origen = os.path.join(RUTA_RESULTADOS_TABLAS, archivo)
                destino = os.path.join(dest_tablas, archivo)
                shutil.move(origen, destino)
                movidos += 1

    # Mover figuras nuevas
    if os.path.exists(RUTA_RESULTADOS_FIGURAS):
        for archivo in os.listdir(RUTA_RESULTADOS_FIGURAS):
            if archivo not in antes_figuras:
                origen = os.path.join(RUTA_RESULTADOS_FIGURAS, archivo)
                destino = os.path.join(dest_figuras, archivo)
                shutil.move(origen, destino)
                movidos += 1

    if movidos > 0:
        print(f"[RESULTADOS] {movidos} archivos movidos a {genome_basename}/")


# =============================================================================
# SERVIDOR HTTP
# =============================================================================

class GenomeHubHandler(http.server.SimpleHTTPRequestHandler):
    """Handler que sirve archivos estáticos y endpoints API."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=DIRECTORIO_BASE, **kwargs)

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        query = urllib.parse.parse_qs(parsed.query)

        # --- ENDPOINTS API ---
        if path == "/api/search":
            q = query.get("query", [""])[0]
            limit = int(query.get("limit", ["20"])[0])
            if not q or len(q) < 3:
                self._json_response({"success": False, "error": "Query debe tener al menos 3 caracteres"}, 400)
                return
            print(f"[BUSCAR] Buscando en NCBI: '{q}' (limit={limit})")
            resultado = buscar_ncbi(q, limit)
            print(f"[BUSCAR] {resultado.get('count', 0)} resultados encontrados")
            self._json_response(resultado)
            return

        if path == "/api/genomes":
            resultado = listar_genomas()
            self._json_response(resultado)
            return

        if path == "/api/results":
            tipo = query.get("type", ["tablas"])[0]
            genome = query.get("genome", [None])[0]
            resultado = listar_resultados(tipo, genome)
            self._json_response(resultado)
            return

        if path == "/api/results_genomes":
            resultado = listar_genomas_con_resultados()
            self._json_response(resultado)
            return

        if path == "/api/generar_informe":
            genome = query.get("genome", [""])[0]
            if not genome:
                self._json_response({"success": False, "error": "Falta el genoma"}, 400)
                return
            # Ejecutar script de generacion de informe
            resultado = ejecutar_analisis("generar_informe", genome_basename=genome)
            if resultado.get("success") and resultado.get("return_code") == 0:
                # Verificar que el PDF se genero
                pdf_path = os.path.join(RUTA_RESULTADOS, genome, f"informe_{genome}.pdf")
                if os.path.exists(pdf_path):
                    # Servir el PDF como descarga
                    self.send_response(200)
                    self.send_header("Content-Type", "application/pdf")
                    self.send_header("Content-Disposition", f'attachment; filename="informe_{genome}.pdf"')
                    self.send_header("Content-Length", str(os.path.getsize(pdf_path)))
                    self.end_headers()
                    with open(pdf_path, "rb") as f:
                        self.wfile.write(f.read())
                    return
                else:
                    self._json_response({"success": False, "error": "El PDF no se genero correctamente", "output": resultado.get("output", "")})
            else:
                self._json_response({"success": False, "error": resultado.get("error", "Error al generar informe"), "output": resultado.get("output", "")})
            return

        if path == "/api/buscar_secuencia":
            genome = query.get("genome", [""])[0]
            secuencia = query.get("secuencia", [""])[0].upper().strip()
            if not genome or not secuencia:
                self._json_response({"success": False, "error": "Faltan genome y secuencia"}, 400)
                return
            if len(secuencia) < 4:
                self._json_response({"success": False, "error": "La secuencia debe tener al menos 4 nucleotidos"}, 400)
                return
            resultado = buscar_secuencia_en_genoma(genome, secuencia)
            self._json_response(resultado)
            return

        if path == "/api/result_data":
            genome = query.get("genome", [""])[0]
            filename = query.get("file", [""])[0]
            if not genome or not filename:
                self._json_response({"success": False, "error": "Faltan genome y file"}, 400)
                return
            if ".." in filename or "/" in filename or "\\" in filename:
                self._json_response({"success": False, "error": "Nombre inválido"}, 400)
                return
            resultado = leer_resultado(genome, filename)
            self._json_response(resultado)
            return

        # Ruta raíz -> app.html
        if path == "/" or path == "":
            self.path = "/app.html"

        # Servir archivos estáticos
        super().do_GET()

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path == "/api/download":
            content_length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(content_length).decode('utf-8')

            try:
                data = json.loads(body)
            except json.JSONDecodeError:
                self._json_response({"success": False, "error": "JSON inválido"}, 400)
                return

            accession_id = data.get("accession_id", "")
            organism_name = data.get("organism_name", "")

            if not accession_id or not organism_name:
                self._json_response({"success": False, "error": "Faltan accession_id y organism_name"}, 400)
                return

            print(f"[DESCARGAR] Descargando {organism_name} ({accession_id})...")
            resultado = descargar_genoma(accession_id, organism_name)

            if resultado["success"]:
                print(f"[DESCARGAR] Descarga completada")
            else:
                print(f"[DESCARGAR] Error: {resultado['error']}")

            self._json_response(resultado)
            return

        if path == "/api/run_analysis":
            content_length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(content_length).decode('utf-8')

            try:
                data = json.loads(body)
            except json.JSONDecodeError:
                self._json_response({"success": False, "error": "JSON inválido"}, 400)
                return

            script = data.get("script", "")
            organism = data.get("organism")
            genome_basename = data.get("genome_basename")
            genome_basename_2 = data.get("genome_basename_2")

            if not script:
                self._json_response({"success": False, "error": "Falta el parámetro 'script'"}, 400)
                return

            print(f"[ANÁLISIS] Ejecutando {script} (genoma: {genome_basename})...")
            resultado = ejecutar_analisis(script, organism, genome_basename, genome_basename_2)

            if resultado.get("success"):
                print(f"[ANÁLISIS] Completado en {resultado.get('execution_time', '?')}s")
            else:
                print(f"[ANÁLISIS] Error: {resultado.get('error', '?')}")

            self._json_response(resultado)
            return

        if path == "/api/delete_genome":
            content_length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(content_length).decode('utf-8')
            try:
                data = json.loads(body)
            except json.JSONDecodeError:
                self._json_response({"success": False, "error": "JSON inválido"}, 400)
                return

            basename = data.get("basename", "")
            if not basename or ".." in basename or "/" in basename or "\\" in basename:
                self._json_response({"success": False, "error": "Nombre inválido"}, 400)
                return

            resultado = eliminar_genoma(basename)
            self._json_response(resultado)
            return

        if path == "/api/delete_result":
            content_length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(content_length).decode('utf-8')
            try:
                data = json.loads(body)
            except json.JSONDecodeError:
                self._json_response({"success": False, "error": "JSON inválido"}, 400)
                return

            filename = data.get("filename", "")
            tipo = data.get("type", "tablas")
            genome = data.get("genome", "")
            if not filename or ".." in filename or "/" in filename or "\\" in filename:
                self._json_response({"success": False, "error": "Nombre inválido"}, 400)
                return

            resultado = eliminar_resultado(filename, tipo, genome)
            self._json_response(resultado)
            return

        if path == "/api/delete_results_all":
            content_length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(content_length).decode('utf-8')
            try:
                data = json.loads(body)
            except json.JSONDecodeError:
                data = {}

            genome = data.get("genome", "")
            resultado = eliminar_todos_resultados_genoma(genome)
            self._json_response(resultado)
            return

        if path == "/api/chat":
            content_length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(content_length).decode('utf-8')
            try:
                data = json.loads(body)
            except json.JSONDecodeError:
                self._json_response({"success": False, "error": "JSON invalido"}, 400)
                return

            message = data.get("message", "")
            context = data.get("context", "")
            system_prompt = data.get("system", "")

            if not message:
                self._json_response({"success": False, "error": "Falta el mensaje"}, 400)
                return

            prompt_completo = message
            if system_prompt:
                prompt_completo = f"{system_prompt}\n\n{prompt_completo}"

            try:
                from backend.gemini_client import consultar_gemini
                respuesta, error = consultar_gemini(prompt_completo, context)
                if respuesta:
                    self._json_response({"success": True, "response": respuesta})
                else:
                    self._json_response({"success": False, "error": error or "Sin respuesta"})
            except Exception as e:
                self._json_response({"success": False, "error": str(e)})
            return

        self._json_response({"error": "Endpoint no encontrado"}, 404)

    def _json_response(self, data, status=200):
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(json.dumps(data, ensure_ascii=False, indent=2).encode("utf-8"))

    def do_OPTIONS(self):
        self.send_response(200)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def log_message(self, format, *args):
        msg = format % args
        if "/api/" in msg or "GET / " in msg:
            print(f"[HTTP] {msg}")


# =============================================================================
# FUNCIONES DE ELIMINACIÓN
# =============================================================================

def eliminar_genoma(basename):
    """Elimina un genoma descargado y sus archivos asociados."""
    archivos = [
        os.path.join(RUTA_DATOS_CRUDO, f"{basename}.gb"),
        os.path.join(RUTA_DATOS_CRUDO, f"{basename}.fasta"),
        os.path.join(RUTA_DATOS_CRUDO, f"metadata_{basename}.json"),
    ]
    eliminados = 0
    for ruta in archivos:
        if os.path.exists(ruta):
            os.remove(ruta)
            eliminados += 1

    # También eliminar resultados de este genoma
    ruta_resultados_genoma = os.path.join(RUTA_RESULTADOS, basename)
    if os.path.exists(ruta_resultados_genoma):
        shutil.rmtree(ruta_resultados_genoma)
        print(f"[ELIMINAR] Resultados de '{basename}' eliminados")

    if eliminados == 0:
        return {"success": False, "error": "No se encontraron archivos para este genoma"}
    print(f"[ELIMINAR] Genoma '{basename}' eliminado ({eliminados} archivos)")
    return {"success": True, "deleted": eliminados}


def eliminar_resultado(filename, tipo="tablas", genome=""):
    """Elimina un archivo de resultado específico."""
    if genome:
        ruta = os.path.join(RUTA_RESULTADOS, genome, tipo, filename)
    else:
        ruta = os.path.join(RUTA_RESULTADOS, tipo, filename)
    if not os.path.exists(ruta):
        return {"success": False, "error": "Archivo no encontrado"}
    os.remove(ruta)
    print(f"[ELIMINAR] Resultado '{filename}' eliminado")
    return {"success": True}


def eliminar_todos_resultados_genoma(genome):
    """Elimina todos los resultados de un genoma."""
    if not genome:
        return {"success": False, "error": "Falta el genoma"}

    ruta = os.path.join(RUTA_RESULTADOS, genome)
    if not os.path.exists(ruta):
        return {"success": False, "error": "No hay resultados para este genoma"}

    total = 0
    for root, dirs, files in os.walk(ruta):
        for f in files:
            os.remove(os.path.join(root, f))
            total += 1

    print(f"[ELIMINAR] {total} resultados de '{genome}' eliminados")
    return {"success": True, "deleted": total}


# =============================================================================
# FUNCIONES DE RESULTADOS
# =============================================================================

def listar_genomas_con_resultados():
    """Lista los genomas que tienen carpetas de resultados."""
    genomas = []
    if not os.path.exists(RUTA_RESULTADOS):
        return {"success": True, "genomes": []}

    for nombre in sorted(os.listdir(RUTA_RESULTADOS)):
        ruta = os.path.join(RUTA_RESULTADOS, nombre)
        # Solo carpetas que NO sean 'tablas' o 'figuras' (legacy)
        if os.path.isdir(ruta) and nombre not in ("tablas", "figuras"):
            # Contar archivos
            num_tablas = 0
            num_figuras = 0
            ruta_t = os.path.join(ruta, "tablas")
            ruta_f = os.path.join(ruta, "figuras")
            if os.path.exists(ruta_t):
                num_tablas = len([f for f in os.listdir(ruta_t) if os.path.isfile(os.path.join(ruta_t, f))])
            if os.path.exists(ruta_f):
                num_figuras = len([f for f in os.listdir(ruta_f) if os.path.isfile(os.path.join(ruta_f, f))])

            if num_tablas > 0 or num_figuras > 0:
                genomas.append({
                    "basename": nombre,
                    "label": nombre.replace("_", " ").title(),
                    "tablas": num_tablas,
                    "figuras": num_figuras
                })

    return {"success": True, "genomes": genomas}


def listar_resultados(tipo="tablas", genome=None):
    """Lista archivos de resultados, filtrado por genoma."""
    resultados = []

    if genome:
        directorio = os.path.join(RUTA_RESULTADOS, genome, tipo)
        path_prefix = f"backend/resultados/{genome}/{tipo}"
    else:
        directorio = os.path.join(RUTA_RESULTADOS, tipo)
        path_prefix = f"backend/resultados/{tipo}"

    if not os.path.exists(directorio):
        return {"success": True, "type": tipo, "genome": genome, "count": 0, "results": []}

    extensiones = {"tablas": [".json", ".csv"], "figuras": [".png"]}
    exts = extensiones.get(tipo, [".json"])

    for archivo in sorted(os.listdir(directorio)):
        if not any(archivo.endswith(e) for e in exts):
            continue

        ruta = os.path.join(directorio, archivo)
        resultados.append({
            "filename": archivo,
            "path": f"{path_prefix}/{archivo}",
            "extension": os.path.splitext(archivo)[1][1:],
            "size_kb": round(os.path.getsize(ruta) / 1024, 2),
            "modified": datetime.fromtimestamp(os.path.getmtime(ruta)).strftime("%Y-%m-%d %H:%M:%S")
        })

    return {"success": True, "type": tipo, "genome": genome, "count": len(resultados), "results": resultados}


def buscar_secuencia_en_genoma(genome_basename, secuencia):
    """Busca una secuencia de nucleotidos en el genoma y devuelve posiciones."""
    archivo_gb = os.path.join(RUTA_DATOS_CRUDO, f"{genome_basename}.gb")
    if not os.path.exists(archivo_gb):
        return {"success": False, "error": "Genoma no encontrado"}

    try:
        from Bio.Seq import Seq
        registro = SeqIO.read(archivo_gb, "genbank")
        seq_str = str(registro.seq).upper()
        secuencia = secuencia.upper()

        # Construir indice de genes
        genes_idx = []
        for feature in registro.features:
            if feature.type == "CDS":
                genes_idx.append({
                    "start": int(feature.location.start),
                    "end": int(feature.location.end),
                    "strand": "+" if feature.location.strand == 1 else "-",
                    "gene": feature.qualifiers.get("gene", [""])[0],
                    "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
                    "product": feature.qualifiers.get("product", [""])[0]
                })

        def find_gene_context(pos_start, pos_end):
            for g in genes_idx:
                if g["start"] <= pos_start < g["end"] or g["start"] < pos_end <= g["end"]:
                    return {
                        "nombre": g["gene"] or g["locus_tag"],
                        "producto": g["product"],
                        "hebra": g["strand"],
                        "inicio": g["start"],
                        "fin": g["end"]
                    }
            return {"nombre": "Intergenico", "producto": "Fuera de CDS", "hebra": "N/A"}

        matches = []

        # Buscar en hebra forward
        start = 0
        while True:
            pos = seq_str.find(secuencia, start)
            if pos == -1:
                break
            context_gene = find_gene_context(pos, pos + len(secuencia))
            context_seq = seq_str[max(0, pos - 20):pos + len(secuencia) + 20]
            matches.append({
                "posicion": pos + 1,
                "fin": pos + len(secuencia),
                "hebra": "+",
                "gen": context_gene,
                "contexto": context_seq
            })
            start = pos + 1
            if len(matches) > 500:
                break

        # Buscar reverse complement
        rev_comp = str(Seq(secuencia).reverse_complement())
        start = 0
        while True:
            pos = seq_str.find(rev_comp, start)
            if pos == -1:
                break
            context_gene = find_gene_context(pos, pos + len(rev_comp))
            context_seq = seq_str[max(0, pos - 20):pos + len(rev_comp) + 20]
            matches.append({
                "posicion": pos + 1,
                "fin": pos + len(rev_comp),
                "hebra": "-",
                "gen": context_gene,
                "contexto": context_seq
            })
            start = pos + 1
            if len(matches) > 500:
                break

        return {
            "success": True,
            "secuencia_buscada": secuencia,
            "longitud_genoma": len(seq_str),
            "total_matches": len(matches),
            "matches": matches[:200]
        }
    except Exception as e:
        return {"success": False, "error": str(e)}


def leer_resultado(genome, filename):
    """Lee el contenido de un archivo de resultado y lo retorna como JSON."""
    import csv as csv_mod

    ruta = os.path.join(RUTA_RESULTADOS, genome, "tablas", filename)
    if not os.path.exists(ruta):
        return {"success": False, "error": "Archivo no encontrado"}

    try:
        if filename.endswith(".json"):
            with open(ruta, "r", encoding="utf-8") as f:
                data = json.load(f)
            return {"success": True, "type": "json", "data": data}

        elif filename.endswith(".csv"):
            rows = []
            with open(ruta, "r", encoding="utf-8") as f:
                reader = csv_mod.DictReader(f)
                headers = reader.fieldnames or []
                for row in reader:
                    rows.append(row)
            return {"success": True, "type": "csv", "headers": headers, "data": rows}

        else:
            return {"success": False, "error": "Tipo de archivo no soportado"}

    except Exception as e:
        return {"success": False, "error": f"Error al leer: {str(e)}"}


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("  GenomeHub - Servidor Local")
    print("=" * 60)
    print(f"  Puerto:    {PUERTO}")
    print(f"  Directorio: {DIRECTORIO_BASE}")
    print(f"  Genomas:   datos/crudo/")
    print(f"  Email NCBI: {CORREO_ELECTRONICO}")
    print()
    print(f"  Abre en tu navegador:")
    print(f"  >>> http://localhost:{PUERTO}")
    print()
    print("  Presiona Ctrl+C para detener")
    print("=" * 60)

    import socket as _socket

    class ReusableHTTPServer(http.server.HTTPServer):
        allow_reuse_address = True
        def server_bind(self):
            self.socket.setsockopt(_socket.SOL_SOCKET, _socket.SO_REUSEADDR, 1)
            super().server_bind()

    try:
        server = ReusableHTTPServer(("", PUERTO), GenomeHubHandler)
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n[INFO] Servidor detenido")
        server.server_close()

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
        # Intentar varios nombres posibles de metadata
        ruta_meta = None
        posibles_meta = [
            os.path.join(RUTA_DATOS_CRUDO, f"metadata_{basename}.json"),
            os.path.join(DIRECTORIO_BASE, "backend", "crudo", f"metadata_{basename}.json"),
            # Nombres cortos (metadata_ecoli.json, metadata_salmonella.json)
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


def ejecutar_analisis(script_name, organism=None):
    """Ejecuta un script de análisis Python."""
    scripts_permitidos = {
        "analisis_genes": {"file": "analisis_genes.py", "requires_organism": True, "timeout": 120},
        "analisis_codones": {"file": "analisis_codones.py", "requires_organism": False, "timeout": 60},
        "comparar_genomas": {"file": "comparar_genomas.py", "requires_organism": False, "timeout": 120},
        "visualizaciones": {"file": "visualizaciones.py", "requires_organism": False, "timeout": 180},
        "analisis_distancias_intergenicas": {"file": "analisis_distancias_intergenicas.py", "requires_organism": False, "timeout": 60},
    }

    if script_name not in scripts_permitidos:
        return {"success": False, "error": f"Script no permitido. Válidos: {', '.join(scripts_permitidos.keys())}"}

    config = scripts_permitidos[script_name]
    script_path = os.path.join(DIRECTORIO_BASE, "backend", "scripts", config["file"])

    if not os.path.exists(script_path):
        return {"success": False, "error": f"Script no encontrado: {config['file']}"}

    # Construir comando
    cmd = [sys.executable, script_path]

    # Mapear organism a número si es necesario
    if config["requires_organism"]:
        organism_map = {"ecoli_k12": "1", "salmonella_lt2": "2"}
        if not organism or organism not in organism_map:
            return {"success": False, "error": "Se requiere organism: ecoli_k12 o salmonella_lt2"}
        cmd.append(organism_map[organism])

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
            resultado = listar_resultados(tipo)
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
            # Leer body JSON
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
                print(f"[DESCARGAR] ✓ Descarga completada -> datos/crudo/")
            else:
                print(f"[DESCARGAR] ✕ Error: {resultado['error']}")

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

            if not script:
                self._json_response({"success": False, "error": "Falta el parámetro 'script'"}, 400)
                return

            print(f"[ANÁLISIS] Ejecutando {script}...")
            resultado = ejecutar_analisis(script, organism)

            if resultado.get("success"):
                print(f"[ANÁLISIS] ✓ Completado en {resultado.get('execution_time', '?')}s")
            else:
                print(f"[ANÁLISIS] ✕ Error: {resultado.get('error', '?')}")

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
            if not filename or ".." in filename or "/" in filename or "\\" in filename:
                self._json_response({"success": False, "error": "Nombre inválido"}, 400)
                return

            resultado = eliminar_resultado(filename, tipo)
            self._json_response(resultado)
            return

        if path == "/api/delete_results_all":
            content_length = int(self.headers.get('Content-Length', 0))
            body = self.rfile.read(content_length).decode('utf-8')
            try:
                data = json.loads(body)
            except json.JSONDecodeError:
                data = {}

            tipo = data.get("type")
            resultado = eliminar_todos_resultados(tipo)
            self._json_response(resultado)
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
        # Solo loguear requests importantes, no cada archivo estático
        msg = format % args
        if "/api/" in msg or "GET / " in msg:
            print(f"[HTTP] {msg}")


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
    if eliminados == 0:
        return {"success": False, "error": "No se encontraron archivos para este genoma"}
    print(f"[ELIMINAR] Genoma '{basename}' eliminado ({eliminados} archivos)")
    return {"success": True, "deleted": eliminados}


def eliminar_resultado(filename, tipo="tablas"):
    """Elimina un archivo de resultado específico."""
    ruta = os.path.join(RUTA_RESULTADOS, tipo, filename)
    if not os.path.exists(ruta):
        return {"success": False, "error": "Archivo no encontrado"}
    os.remove(ruta)
    print(f"[ELIMINAR] Resultado '{filename}' eliminado")
    return {"success": True}


def eliminar_todos_resultados(tipo=None):
    """Elimina todos los resultados, opcionalmente filtrado por tipo."""
    tipos = [tipo] if tipo else ["tablas", "figuras"]
    total = 0
    for t in tipos:
        directorio = os.path.join(RUTA_RESULTADOS, t)
        if not os.path.exists(directorio):
            continue
        for archivo in os.listdir(directorio):
            ruta = os.path.join(directorio, archivo)
            if os.path.isfile(ruta):
                os.remove(ruta)
                total += 1
    print(f"[ELIMINAR] {total} resultados eliminados")
    return {"success": True, "deleted": total}


def listar_resultados(tipo="tablas"):
    """Lista archivos de resultados."""
    resultados = []
    directorio = os.path.join(RUTA_RESULTADOS, tipo)

    if not os.path.exists(directorio):
        return {"success": True, "type": tipo, "count": 0, "results": []}

    extensiones = {"tablas": [".json", ".csv"], "figuras": [".png"]}
    exts = extensiones.get(tipo, [".json"])

    for archivo in sorted(os.listdir(directorio)):
        if not any(archivo.endswith(e) for e in exts):
            continue

        ruta = os.path.join(directorio, archivo)
        resultados.append({
            "filename": archivo,
            "path": f"backend/resultados/{tipo}/{archivo}",
            "extension": os.path.splitext(archivo)[1][1:],
            "size_kb": round(os.path.getsize(ruta) / 1024, 2),
            "modified": datetime.fromtimestamp(os.path.getmtime(ruta)).strftime("%Y-%m-%d %H:%M:%S")
        })

    return {"success": True, "type": tipo, "count": len(resultados), "results": resultados}


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

    try:
        server = http.server.HTTPServer(("", PUERTO), GenomeHubHandler)
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n[INFO] Servidor detenido")
        server.server_close()

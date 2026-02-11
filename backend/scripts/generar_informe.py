#!/usr/bin/env python3
"""
GenomeHub - Generador de Informe PDF en formato IEEE

Recopila todos los analisis de un genoma, genera interpretaciones con IA,
y produce un PDF academico de 20-40 paginas.
"""

import os
import sys
import json
import glob
from datetime import datetime

# ReportLab imports
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch, cm
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle,
    PageBreak, KeepTogether, HRFlowable
)
from reportlab.platypus.tableofcontents import TableOfContents
from reportlab.lib import colors

# =============================================================================
# CONFIGURACION
# =============================================================================

DIRECTORIO_BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # backend/
DIRECTORIO_PROYECTO = os.path.dirname(DIRECTORIO_BASE)  # raiz
RUTA_RESULTADOS = os.path.join(DIRECTORIO_BASE, "resultados")
RUTA_DATOS_CRUDO = os.path.join(DIRECTORIO_PROYECTO, "datos", "crudo")

if len(sys.argv) < 2:
    print("[ERROR] Uso: python generar_informe.py <genome_basename>")
    sys.exit(1)

GENOME_BASENAME = sys.argv[1]


# =============================================================================
# RECOPILAR DATOS
# =============================================================================

def recopilar_datos(genome):
    """Carga todos los JSON de analisis y lista figuras PNG."""
    dir_tablas = os.path.join(RUTA_RESULTADOS, genome, "tablas")
    dir_figuras = os.path.join(RUTA_RESULTADOS, genome, "figuras")

    datos = {}
    figuras = []

    # Cargar JSONs
    if os.path.exists(dir_tablas):
        for f in sorted(os.listdir(dir_tablas)):
            if f.endswith(".json"):
                ruta = os.path.join(dir_tablas, f)
                try:
                    with open(ruta, "r", encoding="utf-8") as fh:
                        datos[f.replace(".json", "")] = json.load(fh)
                    print(f"  [OK] Cargado: {f}")
                except Exception as e:
                    print(f"  [WARN] Error en {f}: {e}")

    # Listar figuras
    if os.path.exists(dir_figuras):
        for f in sorted(os.listdir(dir_figuras)):
            if f.endswith(".png"):
                figuras.append(os.path.join(dir_figuras, f))
                print(f"  [OK] Figura: {f}")

    return datos, figuras


# =============================================================================
# GENERAR INTERPRETACIONES CON IA
# =============================================================================

def generar_interpretaciones_ia(datos, genome):
    """Genera interpretaciones de cada seccion via Gemini."""
    interpretaciones = {}

    try:
        sys.path.insert(0, DIRECTORIO_PROYECTO)
        from backend.gemini_client import consultar_gemini
    except Exception as e:
        print(f"  [WARN] No se pudo cargar gemini_client: {e}")
        return interpretaciones

    nombre_org = genome.replace("_", " ").title()

    # Buscar nombre real en los datos
    for key, data in datos.items():
        if isinstance(data, dict) and data.get("organismo"):
            nombre_org = data["organismo"]
            break

    secciones = {
        "abstract": (
            f"Escribe, únicamente, el abstract en español (150-200 palabras) para un informe científico en formato IEEE sobre {nombre_org}. "
            "Escribe como un investigador en genómica/biólogo molecular: tono técnico y profesional, sin frases introductorias ni comentarios meta (no incluyas 'Claro, aquí tienes' ni 'A continuación'). "
            "Incluye objetivo, método breve (mencionar análisis bioinformático con BioPython), resultados cuantitativos clave (usar los valores proporcionados en el contexto de datos cuando existan) y una oración de conclusión. "
            "Entrega el texto como un único párrafo listo para pegar en el PDF."
        ),

        "introduccion": (
            f"Redacta una introducción académica en español (400-500 palabras) sobre {nombre_org} escrita desde la voz de un biólogo molecular / ingeniero en genómica. "
            "No incluyas avisos ni explicaciones meta; presenta solo el texto final. Describe la importancia biológica del organismo, el contexto científico, el objetivo del estudio, y una visión general del enfoque metodológico (mencionar BioPython y análisis bioinformático). "
            "Evita frases informales; mantén un registro académico apropiado para un informe técnico."
        ),

        "discusion": (
            f"Escribe una discusión técnica en español (400-500 palabras) sobre los resultados del análisis genómico de {nombre_org}. "
            "Actúa como un experto en genómica: interpreta los valores numéricos proporcionados, compara con literatura cuando sea relevante, discute implicaciones biológicas y limitaciones metodológicas, y propone líneas de trabajo futuras. "
            "No abras con prefacios ni justificaciones; entrega solamente el texto de la discusión con tono académico y citas entre corchetes si lo consideras necesario."
        ),

        "conclusiones": (
            f"Genera 5 a 7 conclusiones concisas en español sobre {nombre_org}, cada una en una línea separada. "
            "Escribe desde la perspectiva de un investigador (biólogo/ingeniero), basadas en los datos suministrados; evita frases generales o prefacios. "
            "Cada conclusión debe ser una oración clara que pueda incluir valores numéricos concretos si están disponibles."
        )
    }

    # Agregar contexto de datos reales si hay
    contexto_datos = ""
    for key, data in datos.items():
        if "analisis_genes" in key and isinstance(data, dict):
            stats = data.get("estadisticas_generales", {})
            if stats:
                contexto_datos += f"\nDatos reales: {stats.get('total_genes', 'N/A')} genes CDS, "
                contexto_datos += f"GC {stats.get('contenido_gc_cds', {}).get('promedio', 'N/A')}%, "
                contexto_datos += f"densidad genica {stats.get('densidad_genica_porcentaje', 'N/A')}%"
            gc_pos = data.get("gc_por_posicion_codon", {})
            if gc_pos:
                contexto_datos += f"\nGC por posicion codon: GC1={gc_pos.get('gc1_porcentaje', 'N/A')}%, GC2={gc_pos.get('gc2_porcentaje', 'N/A')}%, GC3={gc_pos.get('gc3_porcentaje', 'N/A')}%"
            esenciales = data.get("genes_esenciales", {})
            if esenciales:
                contexto_datos += f"\nGenes esenciales: {esenciales.get('total_encontrados', 0)}/{esenciales.get('total_referencia', 0)} ({esenciales.get('porcentaje_encontrados', 0)}% cobertura Keio)"
        if "analisis_codones" in key and isinstance(data, dict):
            rscu = data.get("rscu", {})
            if rscu:
                contexto_datos += f"\nRSCU: {rscu.get('total_preferidos', 0)} codones preferidos, {rscu.get('total_evitados', 0)} evitados, {rscu.get('total_raros', 0)} raros"
            nc = data.get("numero_efectivo_codones", {})
            if nc:
                contexto_datos += f"\nNc (numero efectivo de codones): {nc.get('nc', nc.get('valor', 'N/A'))}"

    for seccion, prompt in secciones.items():
        print(f"  [IA] Generando {seccion}...")
        prompt_completo = prompt
        if contexto_datos:
            prompt_completo += f"\n\nContexto con datos reales del analisis:{contexto_datos}"

        respuesta, error = consultar_gemini(prompt_completo)
        if respuesta:
            interpretaciones[seccion] = respuesta
            print(f"  [OK] {seccion} generado ({len(respuesta)} chars)")
        else:
            print(f"  [WARN] No se pudo generar {seccion}: {error}")
            interpretaciones[seccion] = f"[Seccion {seccion} - contenido no disponible]"

    return interpretaciones


def generar_preguntas_ia(datos, genome):
    """Genera preguntas y respuestas para los 2 JSONs mas importantes usando Gemini.
    Devuelve un diccionario {clave_json: texto_qa}.
    """
    qas = {}
    try:
        sys.path.insert(0, DIRECTORIO_PROYECTO)
        from backend.gemini_client import consultar_gemini
    except Exception as e:
        print(f"  [WARN] No se pudo cargar gemini_client para QA: {e}")
        return qas

    nombre_genoma = genome.replace('_', ' ')

    # Solo los 2 JSONs mas relevantes para no exceder timeout
    claves_prioritarias = []
    for key in datos:
        if "analisis_genes" in key:
            claves_prioritarias.insert(0, key)
        elif "analisis_codones" in key and len(claves_prioritarias) < 2:
            claves_prioritarias.append(key)
        elif "analisis_distancias" in key and len(claves_prioritarias) < 2:
            claves_prioritarias.append(key)
    claves_prioritarias = claves_prioritarias[:2]

    for key in claves_prioritarias:
        data = datos[key]
        prompt = (
            f"Investigador en genomica: genera 5 preguntas y respuestas tecnicas sobre '{key}' "
            f"del genoma {nombre_genoma}. Tono academico, sin prefacios. Basate en los datos.\n\n"
            f"Datos:\n{json.dumps(data, indent=2)[:4000]}"
        )

        print(f"  [IA-QA] Generando QA para {key}...")
        respuesta, error = consultar_gemini(prompt)
        if respuesta:
            qas[key] = respuesta
            print(f"  [OK] QA generado para {key} ({len(respuesta)} chars)")
        else:
            print(f"  [WARN] No se pudo generar QA para {key}: {error}")

    return qas


# =============================================================================
# CONSTRUIR PDF
# =============================================================================

def crear_estilos():
    """Crea estilos para el PDF en formato IEEE."""
    styles = getSampleStyleSheet()

    styles.add(ParagraphStyle(
        name='IEEETitle',
        fontName='Times-Bold',
        fontSize=22,
        alignment=TA_CENTER,
        spaceAfter=12,
        textColor=HexColor('#1a1a2e')
    ))

    styles.add(ParagraphStyle(
        name='IEEEAuthor',
        fontName='Times-Roman',
        fontSize=12,
        alignment=TA_CENTER,
        spaceAfter=6,
        textColor=HexColor('#333333')
    ))

    styles.add(ParagraphStyle(
        name='IEEEAbstractTitle',
        fontName='Times-Bold',
        fontSize=10,
        alignment=TA_LEFT,
        spaceAfter=4,
        spaceBefore=12,
        textColor=black
    ))

    styles.add(ParagraphStyle(
        name='IEEEAbstract',
        fontName='Times-Italic',
        fontSize=9,
        alignment=TA_JUSTIFY,
        spaceAfter=12,
        leftIndent=20,
        rightIndent=20,
        textColor=HexColor('#333333')
    ))

    styles.add(ParagraphStyle(
        name='IEEEHeading1',
        fontName='Times-Bold',
        fontSize=12,
        alignment=TA_LEFT,
        spaceBefore=18,
        spaceAfter=8,
        textColor=HexColor('#1a1a2e')
    ))

    styles.add(ParagraphStyle(
        name='IEEEHeading2',
        fontName='Times-Bold',
        fontSize=10,
        alignment=TA_LEFT,
        spaceBefore=12,
        spaceAfter=6,
        textColor=HexColor('#2d2d44')
    ))

    styles.add(ParagraphStyle(
        name='IEEEBody',
        fontName='Times-Roman',
        fontSize=10,
        alignment=TA_JUSTIFY,
        spaceAfter=6,
        leading=14,
        textColor=HexColor('#1a1a1a')
    ))

    styles.add(ParagraphStyle(
        name='IEEECaption',
        fontName='Times-Italic',
        fontSize=8,
        alignment=TA_CENTER,
        spaceAfter=12,
        spaceBefore=4,
        textColor=HexColor('#555555')
    ))

    styles.add(ParagraphStyle(
        name='IEEETableHeader',
        fontName='Times-Bold',
        fontSize=8,
        alignment=TA_CENTER,
        textColor=white
    ))

    styles.add(ParagraphStyle(
        name='IEEETableCell',
        fontName='Times-Roman',
        fontSize=8,
        alignment=TA_LEFT,
        textColor=HexColor('#333333')
    ))

    styles.add(ParagraphStyle(
        name='IEEEReference',
        fontName='Times-Roman',
        fontSize=8,
        alignment=TA_JUSTIFY,
        spaceAfter=4,
        leftIndent=20,
        firstLineIndent=-20,
        textColor=HexColor('#333333')
    ))

    return styles


def fmt_num(valor, decimales=0):
    """Formatea un valor numérico de forma segura con separador de miles.
    Si el valor es 'N/A' o no es un número, retorna 'N/A'."""
    if valor == 'N/A' or valor is None:
        return 'N/A'
    try:
        num_val = float(valor) if isinstance(valor, str) else valor
        if decimales > 0:
            return f"{num_val:,.{decimales}f}"
        else:
            return f"{int(num_val):,}"
    except (ValueError, TypeError):
        return str(valor)


def limpiar_texto_ia(texto):
    """Limpia texto de IA para uso en ReportLab (quita markdown y caracteres no soportados)."""
    if not texto:
        return ""
    texto = texto.replace("**", "").replace("*", "")
    texto = texto.replace("```", "").replace("`", "")
    texto = texto.replace("##", "").replace("#", "")
    # Filtrar caracteres no soportados por Times-Roman (emojis, etc)
    texto = "".join(c for c in texto if ord(c) < 0xFFFF and (ord(c) < 0x2600 or ord(c) > 0x27BF)
                    and (ord(c) < 0x1F300 or ord(c) > 0x1FAFF))
    # Escapar caracteres especiales de XML
    texto = texto.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return texto


def construir_pdf(datos, interpretaciones, figuras, genome):
    """Construye el PDF completo en formato IEEE."""
    nombre_org = genome.replace("_", " ").title()
    for key, data in datos.items():
        if isinstance(data, dict) and data.get("organismo"):
            nombre_org = data["organismo"]
            break

    # Directorio de salida
    dir_genome = os.path.join(RUTA_RESULTADOS, genome)
    os.makedirs(dir_genome, exist_ok=True)
    pdf_path = os.path.join(dir_genome, f"informe_{genome}.pdf")

    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=letter,
        topMargin=1 * inch,
        bottomMargin=1 * inch,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch
    )

    styles = crear_estilos()
    elements = []
    fig_counter = [0]
    table_counter = [0]

    def add_heading1(text, num=""):
        prefix = f"{num}. " if num else ""
        elements.append(Paragraph(f"{prefix}{text}", styles['IEEEHeading1']))

    def add_heading2(text, num=""):
        prefix = f"{num} " if num else ""
        elements.append(Paragraph(f"{prefix}{text}", styles['IEEEHeading2']))

    def add_body(text):
        for parrafo in text.split("\n"):
            parrafo = parrafo.strip()
            if parrafo:
                parrafo_limpio = limpiar_texto_ia(parrafo)
                if parrafo_limpio.startswith("- ") or parrafo_limpio.startswith("• "):
                    parrafo_limpio = "  " + parrafo_limpio
                elements.append(Paragraph(parrafo_limpio, styles['IEEEBody']))

    def add_figure(img_path, caption=""):
        fig_counter[0] += 1
        try:
            img = Image(img_path)
            # Escalar para caber en la pagina
            max_w = 6.5 * inch
            max_h = 4 * inch
            w, h = img.drawWidth, img.drawHeight
            ratio = min(max_w / w, max_h / h, 1.0)
            img.drawWidth = w * ratio
            img.drawHeight = h * ratio
            elements.append(Spacer(1, 6))
            elements.append(img)
            cap_text = f"Fig. {fig_counter[0]}. {caption}" if caption else f"Fig. {fig_counter[0]}."
            elements.append(Paragraph(cap_text, styles['IEEECaption']))
        except Exception as e:
            elements.append(Paragraph(f"[Error cargando figura: {e}]", styles['IEEEBody']))

    def add_data_table(headers, rows, caption=""):
        table_counter[0] += 1
        # Header row
        header_paras = [Paragraph(h, styles['IEEETableHeader']) for h in headers]
        data = [header_paras]
        for row in rows[:30]:  # Limitar filas
            data.append([Paragraph(str(v), styles['IEEETableCell']) for v in row])

        col_widths = [max(1.0 * inch, (6.5 * inch) / len(headers))] * len(headers)
        if len(headers) > 5:
            col_widths = None  # Auto

        t = Table(data, colWidths=col_widths, repeatRows=1)
        t.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), HexColor('#1a1a2e')),
            ('TEXTCOLOR', (0, 0), (-1, 0), white),
            ('FONTNAME', (0, 0), (-1, 0), 'Times-Bold'),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
            ('GRID', (0, 0), (-1, -1), 0.5, HexColor('#cccccc')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [white, HexColor('#f5f5f5')]),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
            ('TOPPADDING', (0, 0), (-1, -1), 4),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ]))
        elements.append(Spacer(1, 8))
        elements.append(t)
        if caption:
            elements.append(Paragraph(f"Tabla {table_counter[0]}. {caption}", styles['IEEECaption']))
        elements.append(Spacer(1, 6))

    # =========================================================================
    # PORTADA
    # =========================================================================
    elements.append(Spacer(1, 2 * inch))
    elements.append(Paragraph("Analisis Genomico Integral", styles['IEEETitle']))
    elements.append(Paragraph(nombre_org, styles['IEEETitle']))
    elements.append(Spacer(1, 0.5 * inch))
    elements.append(HRFlowable(width="50%", color=HexColor('#1a1a2e'), thickness=2))
    elements.append(Spacer(1, 0.3 * inch))
    elements.append(Paragraph("Informe de Analisis Bioinformatico", styles['IEEEAuthor']))
    # Solo la fecha (sin prefacio 'Generado por...')
    elements.append(Paragraph(f"{datetime.now().strftime('%d de %B de %Y')}", styles['IEEEAuthor']))
    elements.append(Spacer(1, 0.3 * inch))
    elements.append(Paragraph("Herramientas: BioPython, NCBI GenBank, Gemini AI", styles['IEEEAuthor']))
    elements.append(Paragraph("Formato: IEEE", styles['IEEEAuthor']))
    elements.append(Spacer(1, 0.2 * inch))
    # Bloque institucional requerido por el usuario
    elements.append(Paragraph('Curso: Bioinformatica', styles['IEEEBody']))
    elements.append(Spacer(1, 6))
    elements.append(Paragraph('Carrera: Ing. Informatica y de Sistemas', styles['IEEEBody']))
    elements.append(Paragraph('Universidad: UNSAAC - Cusco, Peru', styles['IEEEBody']))
    elements.append(Paragraph('Facultad de Ingenieria Electrica, Electronica, Informatica y Mecanica', styles['IEEEBody']))
    elements.append(Paragraph('Universidad Nacional San Antonio Abad del Cusco', styles['IEEEBody']))
    elements.append(PageBreak())

    # =========================================================================
    # ABSTRACT
    # =========================================================================
    elements.append(Paragraph("Abstract", styles['IEEEAbstractTitle']))
    abstract_text = limpiar_texto_ia(interpretaciones.get("abstract", "Abstract no disponible."))
    elements.append(Paragraph(abstract_text, styles['IEEEAbstract']))
    elements.append(Spacer(1, 12))
    elements.append(HRFlowable(width="100%", color=HexColor('#cccccc'), thickness=0.5))

    # =========================================================================
    # I. INTRODUCCION
    # =========================================================================
    add_heading1("INTRODUCCION", "I")
    intro_text = interpretaciones.get("introduccion", "Introduccion no disponible.")
    add_body(intro_text)

    # =========================================================================
    # II. MATERIALES Y METODOS
    # =========================================================================
    add_heading1("MATERIALES Y METODOS", "II")

    add_heading2("Obtencion del genoma", "A.")
    add_body(f"El genoma de referencia de {nombre_org} fue obtenido de la base de datos NCBI GenBank en formato GenBank (.gb). El archivo fue descargado mediante la API de Entrez de NCBI utilizando BioPython.")

    add_heading2("Pipeline de analisis bioinformatico", "B.")
    add_body("El analisis se realizo utilizando un pipeline automatizado desarrollado en Python 3 con las siguientes herramientas y librerias:")
    add_body("- BioPython (SeqIO, SeqUtils): Lectura y procesamiento de archivos GenBank, extraccion de features genomicos (CDS, tRNA, rRNA), calculo de contenido GC.")
    add_body("- Pandas y NumPy: Procesamiento estadistico de datos genomicos.")
    add_body("- Matplotlib y Seaborn: Generacion de graficos y visualizaciones.")
    add_body("- Gemini AI (Google): Generacion de interpretaciones biologicas automatizadas.")

    add_heading2("Analisis realizados", "C.")
    add_body("Se realizaron los siguientes analisis sobre el genoma:")
    add_body("1. Analisis de genes: Extraccion de todas las CDS (Coding DNA Sequences), calculo de estadisticas de tamano, contenido GC por gen, distribucion por hebra, y comparacion con valores de literatura.")
    add_body("2. Analisis de codones: Conteo de los 64 codones, frecuencia de uso, codones de inicio y parada, sesgo de codones.")
    add_body("3. Analisis de distancias intergenicas: Calculo de distancias entre genes adyacentes, identificacion de regiones intergenicas grandes (posibles islas genomicas).")
    add_body("4. Analisis de estructura genomica: Composicion del genoma (regiones codificantes vs no codificantes), identificacion de operones putativos, analisis de GC por region.")

    # Agregar tabla de software
    add_data_table(
        ["Software", "Version", "Uso"],
        [
            ["Python", "3.x", "Lenguaje principal"],
            ["BioPython", "1.x", "Analisis de secuencias"],
            ["Matplotlib", "3.x", "Visualizaciones"],
            ["ReportLab", "4.x", "Generacion de PDF"],
            ["Gemini AI", "2.5 Flash Lite", "Interpretaciones IA"],
        ],
        "Herramientas de software utilizadas en el analisis."
    )

    # =========================================================================
    # III. RESULTADOS
    # =========================================================================
    add_heading1("RESULTADOS", "III")

    # --- 3A. Analisis de Genes ---
    genes_key = f"analisis_genes_{genome}"
    if genes_key in datos:
        genes_data = datos[genes_key]
        stats = genes_data.get("estadisticas_generales", {})

        add_heading2("Analisis de Genes Codificantes (CDS)", "A.")
        add_body(f"Se identificaron un total de {fmt_num(stats.get('total_genes', 'N/A'))} secuencias codificantes (CDS) en el genoma de {nombre_org}.")

        if stats:
            gc_cds = stats.get("contenido_gc_cds", {})
            tamano = stats.get("tamano_gen", {})

            add_data_table(
                ["Metrica", "Valor", "Unidad"],
                [
                    ["Total CDS", f"{fmt_num(stats.get('total_genes', 'N/A'))}", "genes"],
                    ["Densidad genica", f"{stats.get('densidad_genica_porcentaje', 'N/A')}", "%"],
                    ["GC promedio en CDS", f"{gc_cds.get('promedio', 'N/A')}", "%"],
                    ["Tamano promedio gen", f"{fmt_num(tamano.get('promedio_pb', 'N/A'), 0)}", "pb"],
                    ["Tamano mediana gen", f"{fmt_num(tamano.get('mediana_pb', 'N/A'), 0)}", "pb"],
                    ["Gen mas grande", f"{fmt_num(tamano.get('maximo_pb', 'N/A'))}", "pb"],
                    ["Gen mas pequeno", f"{fmt_num(tamano.get('minimo_pb', 'N/A'))}", "pb"],
                ],
                "Estadisticas generales de genes codificantes."
            )

        # Genes extremos
        extremos = genes_data.get("genes_extremos", {})
        if extremos:
            mayor = extremos.get("gen_mas_largo", {})
            menor = extremos.get("gen_mas_corto", {})
            add_body(f"El gen mas largo identificado fue {mayor.get('nombre', mayor.get('nombre_gen', 'N/A'))} ({mayor.get('producto', 'N/A')}) con {fmt_num(mayor.get('longitud_pb', 0))} pb, mientras que el mas corto fue {menor.get('nombre', menor.get('nombre_gen', 'N/A'))} ({menor.get('producto', 'N/A')}) con {fmt_num(menor.get('longitud_pb', 0))} pb.")

            # Top 10 mas largos
            top_largos = extremos.get("top_10_mas_largos", [])
            if top_largos:
                add_heading2("Top 10 Genes Mas Largos", "")
                rows_largos = []
                for i, g in enumerate(top_largos, 1):
                    rows_largos.append([
                        str(i),
                        g.get("nombre", "") or g.get("locus_tag", ""),
                        g.get("producto", "Sin anotacion")[:60],
                        f"{fmt_num(g.get('longitud_pb', 0))} pb",
                        str(g.get("num_aminoacidos", 0)) + " aa"
                    ])
                add_data_table(
                    ["#", "Gen", "Proteina", "Tamano", "Aminoacidos"],
                    rows_largos,
                    "Los 10 genes codificantes mas largos del genoma."
                )

            # Top 10 mas cortos
            top_cortos = extremos.get("top_10_mas_cortos", [])
            if top_cortos:
                add_heading2("Top 10 Genes Mas Cortos", "")
                rows_cortos = []
                for i, g in enumerate(top_cortos, 1):
                    rows_cortos.append([
                        str(i),
                        g.get("nombre", "") or g.get("locus_tag", ""),
                        g.get("producto", "Sin anotacion")[:60],
                        f"{fmt_num(g.get('longitud_pb', 0))} pb",
                        str(g.get("num_aminoacidos", 0)) + " aa"
                    ])
                add_data_table(
                    ["#", "Gen", "Proteina", "Tamano", "Aminoacidos"],
                    rows_cortos,
                    "Los 10 genes codificantes mas cortos del genoma."
                )

        # Comparacion con literatura
        lit = genes_data.get("comparacion_literatura", {})
        if lit and len(lit) > 0:
            add_heading2("Comparacion con literatura", "")
            rows_lit = []
            for metrica, info in lit.items():
                if isinstance(info, dict):
                    rows_lit.append([
                        metrica,
                        str(info.get("valor_observado", "N/A")),
                        str(info.get("valor_literatura", "N/A")),
                        str(info.get("diferencia_porcentaje", "N/A")) + "%"
                    ])
            if rows_lit:
                add_data_table(
                    ["Metrica", "Observado", "Literatura", "Diferencia"],
                    rows_lit,
                    "Comparacion de metricas observadas con valores de literatura."
                )

        # GC por Posicion de Codon
        gc_pos = genes_data.get("gc_por_posicion_codon", {})
        if gc_pos:
            add_heading2("GC por Posicion de Codon", "")
            add_body(f"El analisis de contenido GC por posicion del codon mostro: GC1 (1ra posicion) = {gc_pos.get('gc1_porcentaje', 'N/A')}%, GC2 (2da posicion) = {gc_pos.get('gc2_porcentaje', 'N/A')}%, GC3 (3ra posicion) = {gc_pos.get('gc3_porcentaje', 'N/A')}%.")
            if gc_pos.get("interpretacion"):
                add_body(gc_pos["interpretacion"])
            add_data_table(
                ["Posicion", "GC (%)", "Significado"],
                [
                    ["GC1 (1ra posicion)", str(gc_pos.get("gc1_porcentaje", "N/A")) + "%", "Mas conservada, limita aminoacidos"],
                    ["GC2 (2da posicion)", str(gc_pos.get("gc2_porcentaje", "N/A")) + "%", "Afecta directamente al aminoacido"],
                    ["GC3 (3ra posicion)", str(gc_pos.get("gc3_porcentaje", "N/A")) + "%", "Refleja sesgo codonico y presion mutacional"],
                ],
                "Contenido GC por posicion del codon."
            )

        # Genes Esenciales
        esenciales = genes_data.get("genes_esenciales", {})
        if esenciales:
            add_heading2("Genes Esenciales (Keio Collection)", "")
            add_body(f"Se identificaron {esenciales.get('total_encontrados', 0)} de {esenciales.get('total_referencia', 0)} genes esenciales de la coleccion Keio ({esenciales.get('porcentaje_encontrados', 0)}% de cobertura). Fuente: {esenciales.get('fuente', 'Baba et al., 2006')}.")
            cats_esc = esenciales.get("categorias_esenciales", {})
            if cats_esc:
                rows_esc = []
                for cat, count in sorted(cats_esc.items(), key=lambda x: x[1], reverse=True):
                    rows_esc.append([cat, str(count)])
                add_data_table(
                    ["Categoria Funcional", "Genes"],
                    rows_esc[:10],
                    "Categorias funcionales de genes esenciales encontrados."
                )

        # Densidad por Ventana
        densidad = genes_data.get("densidad_por_ventana", {})
        if densidad:
            add_heading2("Densidad Genica por Region", "")
            add_body(f"El analisis de densidad genica por ventana de {fmt_num(densidad.get('ventana_pb', 50000))} pb mostro una densidad promedio de {densidad.get('promedio_densidad', 'N/A')}%, con un maximo de {densidad.get('max_densidad', 'N/A')}% y un minimo de {densidad.get('min_densidad', 'N/A')}%. Se analizaron {densidad.get('total_ventanas', 0)} ventanas a lo largo del genoma.")

    # --- 3B. Analisis de Codones ---
    codones_key = f"analisis_codones_{genome}"
    if codones_key in datos:
        codones_data = datos[codones_key]
        add_heading2("Analisis de Uso de Codones", "B.")

        gc = codones_data.get("contenido_gc", {})
        if gc:
            add_body(f"El contenido GC global del genoma es {gc.get('gc_total_porcentaje', 'N/A')}%. La composicion nucleotidica muestra: A={gc.get('composicion', {}).get('A_porcentaje', 'N/A')}%, T={gc.get('composicion', {}).get('T_porcentaje', 'N/A')}%, G={gc.get('composicion', {}).get('G_porcentaje', 'N/A')}%, C={gc.get('composicion', {}).get('C_porcentaje', 'N/A')}%.")

        inicio = codones_data.get("codones_inicio", {})
        parada = codones_data.get("codones_parada", {})
        if inicio:
            add_body(f"El codon de inicio ATG se utiliza en el {inicio.get('ATG', {}).get('proporcion_porcentaje', 100)}% de los genes ({fmt_num(inicio.get('ATG', {}).get('conteo', 'N/A'))} ocurrencias).")
        if parada:
            codones_stop = []
            for codon in ["TAA", "TAG", "TGA"]:
                if codon in parada:
                    codones_stop.append(f"{codon}: {parada[codon].get('proporcion_porcentaje', 0)}%")
            if codones_stop:
                add_body(f"Distribucion de codones de parada: {', '.join(codones_stop)}.")

        # RSCU
        rscu = codones_data.get("rscu", {})
        if rscu:
            add_heading2("RSCU - Uso Relativo de Codones Sinonimos", "")
            total_pref = rscu.get("total_preferidos", len(rscu.get("codones_preferidos", {})))
            total_evit = rscu.get("total_evitados", len(rscu.get("codones_evitados", {})))
            total_raros = rscu.get("total_raros", len(rscu.get("codones_raros", {})))
            add_body(f"El analisis RSCU (Relative Synonymous Codon Usage) identifico {total_pref} codones preferidos (RSCU > 1.3), {total_evit} codones evitados (RSCU < 0.7) y {total_raros} codones raros (RSCU < 0.5). Un RSCU = 1.0 indica uso sin sesgo entre codones sinonimos.")

            # Tabla de codones preferidos
            preferidos = rscu.get("codones_preferidos", {})
            if preferidos:
                rows_pref = []
                for codon, val in sorted(preferidos.items(), key=lambda x: x[1] if isinstance(x[1], (int, float)) else 0, reverse=True)[:15]:
                    rows_pref.append([codon, str(round(val, 2)) if isinstance(val, (int, float)) else str(val)])
                add_data_table(
                    ["Codon", "RSCU"],
                    rows_pref,
                    "Codones preferidos (RSCU > 1.3) - uso selectivo por eficiencia traduccional."
                )

            # Codones raros
            raros = rscu.get("codones_raros", {})
            if raros:
                raros_list = ", ".join([f"{c}({round(v, 2) if isinstance(v, (int, float)) else v})" for c, v in raros.items()])
                add_body(f"Codones raros (RSCU < 0.5): {raros_list}.")

        # Numero Efectivo de Codones (Nc)
        nc_data = codones_data.get("numero_efectivo_codones", {})
        if nc_data:
            add_heading2("Numero Efectivo de Codones (Nc)", "")
            nc_val = nc_data.get("nc", nc_data.get("valor", 0))
            add_body(f"El Numero Efectivo de Codones (Nc) calculado fue {nc_val} (rango teorico: 20 = sesgo maximo, 61 = sin sesgo). {nc_data.get('interpretacion', '')}. Referencia: {nc_data.get('referencia', 'Wright F. (1990) Gene 87:23-29')}.")

        # Sesgo Codonico por Aminoacido
        sesgo = codones_data.get("sesgo_por_aminoacido", {})
        if sesgo:
            add_heading2("Sesgo Codonico por Aminoacido", "")
            add_body("Para cada aminoacido con codones sinonimos se identifico el codon preferido y el evitado:")
            rows_sesgo = []
            for aa, info in sorted(sesgo.items()):
                if isinstance(info, dict) and info.get("codones") and len(info.get("codones", [])) > 1:
                    rows_sesgo.append([
                        aa,
                        info.get("preferido", "-"),
                        str(info.get("preferido_porcentaje", "-")) + "%",
                        info.get("evitado", "-"),
                        str(info.get("evitado_porcentaje", "-")) + "%"
                    ])
            if rows_sesgo:
                add_data_table(
                    ["Aminoacido", "Preferido", "% Uso", "Evitado", "% Uso"],
                    rows_sesgo,
                    "Sesgo codonico: codon preferido vs evitado por aminoacido."
                )

    # --- 3C. Distancias Intergenicas ---
    dist_key = f"analisis_distancias_{genome}"
    if dist_key in datos:
        dist_data = datos[dist_key]
        add_heading2("Distancias Intergenicas", "C.")

        stats_dist = dist_data.get("estadisticas_generales", {})
        if stats_dist:
            add_body(f"Se analizaron {fmt_num(stats_dist.get('total_pares_adyacentes', 'N/A'))} pares de genes adyacentes. La distancia intergenica promedio fue de {fmt_num(stats_dist.get('distancia_promedio_pb', 0), 1)} pb con una mediana de {fmt_num(stats_dist.get('distancia_mediana_pb', 0), 1)} pb.")
            add_body(f"Se identificaron {fmt_num(stats_dist.get('genes_solapados', 0))} pares de genes solapados ({stats_dist.get('porcentaje_solapados', 0):.1f}% del total) y {stats_dist.get('regiones_grandes_5kb', 0)} regiones intergenicas mayores a 5 kb, que podrian representar islas genomicas o regiones regulatorias.")

        # Top 20 regiones grandes con productos
        top_reg = dist_data.get("top_20_regiones_grandes", [])
        if top_reg:
            rows_reg = []
            for i, r in enumerate(top_reg[:15], 1):
                gen1_prod = r.get("gen1_producto", "")
                gen2_prod = r.get("gen2_producto", "")
                rows_reg.append([
                    str(i),
                    r.get("gen1", ""),
                    (gen1_prod[:40] if gen1_prod else "-"),
                    r.get("gen2", ""),
                    (gen2_prod[:40] if gen2_prod else "-"),
                    f"{fmt_num(r.get('distancia_pb', 0))} pb"
                ])
            add_data_table(
                ["#", "Gen 1", "Proteina Gen 1", "Gen 2", "Proteina Gen 2", "Distancia"],
                rows_reg,
                "Top 15 regiones intergenicas mas grandes (posibles islas genomicas)."
            )

    # --- 3D. Estructura Genomica ---
    est_key = f"analisis_estructura_{genome}"
    if est_key in datos:
        est_data = datos[est_key]
        add_heading2("Estructura Genomica", "D.")

        comp = est_data.get("composicion_genomica", {})
        if comp:
            cds_pct = comp.get("codificante_cds", {}).get("porcentaje", 0)
            inter_pct = comp.get("intergenico", {}).get("porcentaje", 0)
            add_body(f"El genoma presenta una alta densidad codificante tipica de procariotas: {cds_pct}% del genoma corresponde a CDS, mientras que {inter_pct}% son regiones intergenicas.")

        operones = est_data.get("operones_putativos", {})
        if operones:
            add_body(f"Se identificaron {operones.get('total', 0)} operones putativos (grupos de genes adyacentes en la misma hebra con distancia < {operones.get('umbral_distancia_pb', 50)} pb). El {operones.get('porcentaje_genes_en_operones', 0)}% de los genes se encuentran organizados en operones, consistente con la organizacion policistronicas tipica de bacterias.")

        gc_reg = est_data.get("gc_por_region", {})
        if gc_reg:
            add_body(f"El contenido GC difiere entre regiones codificantes ({gc_reg.get('gc_codificante', 0)}%) y no codificantes ({gc_reg.get('gc_no_codificante', 0)}%), con una diferencia de {gc_reg.get('diferencia_gc', 0)} puntos porcentuales.")

    # --- 3E. Analisis de Proteinas ---
    prot_key = f"analisis_proteinas_{genome}"
    if prot_key in datos:
        prot_data = datos[prot_key]
        add_heading2("Analisis de Proteinas", "E.")

        # Estructura Primaria
        primaria = prot_data.get("estructura_primaria", {})
        stats_prot = primaria.get("estadisticas_generales", {})
        if stats_prot:
            add_body(f"Se analizaron {fmt_num(stats_prot.get('total_analizadas', 0))} proteinas del proteoma. La longitud promedio fue de {fmt_num(stats_prot.get('longitud_promedio_aa', 0), 1)} aminoacidos con un peso molecular promedio de {fmt_num(stats_prot.get('peso_molecular_promedio_da', 0), 0)} Da. El punto isoelectrico promedio fue {stats_prot.get('pi_promedio', 0)} y el indice GRAVY promedio fue {stats_prot.get('gravy_promedio', 0)}.")
            add_body(f"Se identificaron {stats_prot.get('proteinas_acidas', 0)} proteinas acidas (pI < 7) y {stats_prot.get('proteinas_basicas', 0)} basicas (pI > 7). {stats_prot.get('proteinas_estables', 0)} proteinas se clasificaron como estables (indice inestabilidad < 40) y {stats_prot.get('proteinas_inestables', 0)} como inestables.")

            add_data_table(
                ["Propiedad", "Valor"],
                [
                    ["Total proteinas", str(stats_prot.get("total_analizadas", 0))],
                    ["Longitud promedio", f"{stats_prot.get('longitud_promedio_aa', 0)} aa"],
                    ["Peso molecular promedio", f"{fmt_num(stats_prot.get('peso_molecular_promedio_da', 0), 0)} Da"],
                    ["pI promedio", str(stats_prot.get("pi_promedio", 0))],
                    ["GRAVY promedio", str(stats_prot.get("gravy_promedio", 0))],
                    ["Proteinas acidas / basicas", f"{stats_prot.get('proteinas_acidas', 0)} / {stats_prot.get('proteinas_basicas', 0)}"],
                    ["Estables / inestables", f"{stats_prot.get('proteinas_estables', 0)} / {stats_prot.get('proteinas_inestables', 0)}"],
                    ["Peptidos senal detectados", str(primaria.get("peptidos_senal", {}).get("total_detectados", 0))],
                ],
                "Estadisticas generales del proteoma."
            )

        # Categorias funcionales
        cats = primaria.get("categorias_funcionales", {})
        if cats:
            rows_cats = []
            for cat_nombre, cat_data in sorted(cats.items(), key=lambda x: x[1].get("total", 0), reverse=True):
                if cat_data.get("total", 0) > 0:
                    rows_cats.append([cat_nombre, str(cat_data["total"]), f"{cat_data.get('porcentaje', 0)}%"])
            if rows_cats:
                add_data_table(
                    ["Categoria", "Total", "Porcentaje"],
                    rows_cats[:10],
                    "Clasificacion funcional de las proteinas del proteoma."
                )

        # Estructura Secundaria
        sec = prot_data.get("estructura_secundaria", {})
        prom_sec = sec.get("promedio_proteoma", {})
        if prom_sec:
            add_body(f"La prediccion de estructura secundaria (metodo Chou-Fasman) mostro que en promedio el proteoma contiene {prom_sec.get('helix', 0)}% helice alfa, {prom_sec.get('sheet', 0)}% lamina beta y {prom_sec.get('turn', 0)}% giros y regiones desordenadas.")

        # Estructura Cuaternaria
        cuat = prot_data.get("estructura_cuaternaria", {})
        if cuat.get("total_complejos", 0) > 0:
            add_body(f"Se detectaron {cuat['total_complejos']} complejos proteicos multi-subunidad con un total de {cuat.get('total_subunidades', 0)} subunidades identificadas por anotacion.")

            top_complejos = cuat.get("complejos_detectados", [])[:10]
            if top_complejos:
                rows_comp = []
                for c in top_complejos:
                    genes = ", ".join([s.get("nombre_gen", "") or s.get("locus_tag", "") for s in c.get("subunidades", [])])
                    rows_comp.append([c.get("nombre_complejo", "")[:45], str(c.get("num_subunidades", 0)), genes[:50]])
                add_data_table(
                    ["Complejo", "Subunidades", "Genes"],
                    rows_comp,
                    "Principales complejos proteicos detectados."
                )

        # Mutaciones patogenicas
        muts = primaria.get("mutaciones_patogenicas", [])
        if muts:
            rows_mut = []
            for m in muts:
                estado = "Presente" if m.get("encontrado") else "No encontrado"
                long_obs = str(m.get("longitud_observada_aa", "-"))
                rows_mut.append([m.get("gen", ""), estado, long_obs, m.get("descripcion", "")[:50]])
            add_data_table(
                ["Gen", "Estado", "Longitud (aa)", "Funcion"],
                rows_mut,
                "Genes criticos para resistencia a antibioticos."
            )

    # --- 3F. Comparacion de Genomas ---
    for key, data in datos.items():
        if "comparacion_" in key and isinstance(data, dict) and data.get("organismos_comparados"):
            add_heading2("Comparacion de Genomas", "F.")

            orgs = data.get("organismos_comparados", {})
            org1 = orgs.get("organismo_1", {}).get("nombre", "Organismo 1")
            org2 = orgs.get("organismo_2", {}).get("nombre", "Organismo 2")
            add_body(f"Se realizo una comparacion genomica entre {org1} y {org2}.")

            metricas = data.get("metricas_generales", {})
            if metricas:
                ec = metricas.get("ecoli", metricas.get(list(metricas.keys())[0], {})) if metricas else {}
                sal = metricas.get("salmonella", metricas.get(list(metricas.keys())[-1], {})) if len(metricas) > 1 else {}

                if ec and sal:
                    add_data_table(
                        ["Metrica", org1[:30], org2[:30]],
                        [
                            ["Longitud (pb)", f"{fmt_num(ec.get('longitud_genoma_pb', 0))}", f"{fmt_num(sal.get('longitud_genoma_pb', 0))}"],
                            ["Total CDS", f"{fmt_num(ec.get('total_genes_cds', 0))}", f"{fmt_num(sal.get('total_genes_cds', 0))}"],
                            ["GC (%)", f"{ec.get('contenido_gc_porcentaje', 0)}", f"{sal.get('contenido_gc_porcentaje', 0)}"],
                            ["Densidad genica (%)", f"{ec.get('densidad_genica_porcentaje', 0)}", f"{sal.get('densidad_genica_porcentaje', 0)}"],
                        ],
                        f"Comparacion de metricas generales entre {org1} y {org2}."
                    )

            # Interpretacion IA de comparacion
            interp_comp = data.get("interpretacion_ia", "")
            if interp_comp:
                add_body(interp_comp)

            break  # Solo una comparacion

    # =========================================================================
    # FIGURAS
    # =========================================================================
    if figuras:
        elements.append(PageBreak())
        add_heading1("FIGURAS", "IV")

        fig_captions = {
            "resumen_general": "Resumen general del analisis genomico",
            "distribucion_tamanos": "Distribucion de tamanos de genes codificantes",
            "distribucion_hebras": "Distribucion de genes por hebra (forward vs reverse)",
            "genes_extremos": "Genes extremos: mas largo y mas corto",
            "genes_vs_cds": "Comparacion genes totales vs CDS",
            "codones_inicio": "Uso del codon de inicio ATG",
            "comparacion_genomas": "Comparacion de metricas entre genomas",
            "comparacion_virulencia": "Comparacion de genes de virulencia",
            "comparacion_distancias": "Comparacion de distancias intergenicas",
            "resumen_comparacion": "Resumen de la comparacion genomica",
        }

        for fig_path in figuras:
            try:
                if os.path.exists(fig_path):
                    fig_name = os.path.splitext(os.path.basename(fig_path))[0]
                    caption = fig_captions.get(fig_name, fig_name.replace("_", " ").title())
                    add_figure(fig_path, caption)
                else:
                    print(f"  [WARN] Figura no encontrada: {fig_path}")
                    elements.append(Paragraph(f"[Figura no encontrada: {os.path.basename(fig_path)}]", styles['IEEEBody']))
            except Exception as e:
                print(f"  [WARN] Error al agregar figura {fig_path}: {e}")
                elements.append(Paragraph(f"[Error al procesar figura: {str(e)}]", styles['IEEEBody']))

    # =========================================================================
    # PREGUNTAS Y RESPUESTAS (IA) POR CADA JSON
    # =========================================================================
    try:
        qas = generar_preguntas_ia(datos, genome)
        if qas:
            elements.append(PageBreak())
            add_heading1("PREGUNTAS Y RESPUESTAS (IA) POR RESULTADO", "III.A")
            for key, texto in qas.items():
                add_heading2(key.replace('_', ' ').upper(), "")
                qa_text = limpiar_texto_ia(texto)
                # separar por lineas y agregar como parrafos
                for par in qa_text.split('\n\n'):
                    if par.strip():
                        elements.append(Paragraph(par.strip(), styles['IEEEBody']))
                        elements.append(Spacer(1, 6))
    except Exception as e:
        print(f"  [WARN] Fallo al generar QA por JSON: {e}")

    # =========================================================================
    # V. DISCUSION
    # =========================================================================
    elements.append(PageBreak())
    add_heading1("DISCUSION", "V")
    discusion_text = interpretaciones.get("discusion", "Discusion no disponible.")
    add_body(discusion_text)

    # =========================================================================
    # VI. CONCLUSIONES
    # =========================================================================
    add_heading1("CONCLUSIONES", "VI")
    conclusiones_text = interpretaciones.get("conclusiones", "Conclusiones no disponibles.")
    add_body(conclusiones_text)

    # =========================================================================
    # REFERENCIAS
    # =========================================================================
    elements.append(PageBreak())
    add_heading1("REFERENCIAS", "")

    refs_list = [
        "[1] Blattner FR et al. (1997). The complete genome sequence of Escherichia coli K-12. Science, 277(5331), 1453-1462.",
        "[2] NCBI GenBank. National Center for Biotechnology Information. https://www.ncbi.nlm.nih.gov/genbank/",
        "[3] Cock PJ et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422-1423.",
        "[4] EcoCyc Database. https://ecocyc.org/",
        "[5] Keseler IM et al. (2021). The EcoCyc Database in 2021. Frontiers in Microbiology, 12, 711077.",
        "[6] Hunter JD (2007). Matplotlib: A 2D graphics environment. Computing in Science &amp; Engineering, 9(3), 90-95.",
    ]

    # Agregar referencias del analisis
    for key, data in datos.items():
        if isinstance(data, dict):
            refs = data.get("referencias_bibliograficas", {})
            if isinstance(refs, dict):
                for ref_key, ref_val in refs.items():
                    if isinstance(ref_val, str) and ref_val not in str(refs_list):
                        refs_list.append(f"[{len(refs_list) + 1}] {ref_val}")

    for ref in refs_list:
        elements.append(Paragraph(ref, styles['IEEEReference']))

    # =========================================================================
    # GENERAR PDF
    # =========================================================================
    print(f"\n  Generando PDF ({len(elements)} elementos)...")
    try:
        doc.build(elements)
        print(f"  [OK] PDF generado: {pdf_path}")
    except Exception as e:
        print(f"  [ERROR] Fallo al generar PDF con ReportLab:")
        print(f"    {type(e).__name__}: {str(e)}")
        import traceback
        traceback.print_exc()
        raise
    return pdf_path


# =============================================================================
# MAIN
# =============================================================================

def main():
    try:
        print("\n" + "=" * 70)
        print(f"GENERACION DE INFORME PDF - FORMATO IEEE")
        print(f"Genoma: {GENOME_BASENAME}")
        print("=" * 70)
        print(f"  Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

        # Paso 1: Recopilar datos
        print("\n[PASO 1/3] Recopilando datos de analisis...")
        datos, figuras = recopilar_datos(GENOME_BASENAME)
        print(f"  Total: {len(datos)} archivos JSON, {len(figuras)} figuras PNG")

        if not datos:
            print("\n[ERROR] No se encontraron datos de analisis. Ejecuta los analisis primero.")
            sys.exit(1)

        # Paso 2: Generar interpretaciones IA
        print("\n[PASO 2/3] Generando interpretaciones con IA...")
        interpretaciones = generar_interpretaciones_ia(datos, GENOME_BASENAME)

        # Paso 3: Construir PDF
        print("\n[PASO 3/3] Construyendo PDF...")
        pdf_path = construir_pdf(datos, interpretaciones, figuras, GENOME_BASENAME)

        file_size = os.path.getsize(pdf_path)
        print(f"\n" + "=" * 70)
        print(f"[OK] Informe PDF generado exitosamente")
        print(f"     Archivo: informe_{GENOME_BASENAME}.pdf")
        print(f"     Tamano: {file_size / 1024:.1f} KB")
        print("=" * 70 + "\n")
    except Exception as e:
        print(f"\n[ERROR] Excepcion durante la generacion del informe:")
        print(f"  {type(e).__name__}: {str(e)}")
        import traceback
        print("\n[TRACEBACK]")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

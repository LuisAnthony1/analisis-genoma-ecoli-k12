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
        "abstract": f"Genera un abstract academico (150-200 palabras) en espanol para un informe de analisis genomico de {nombre_org}. Incluye objetivo, metodos (analisis bioinformatico con BioPython), principales hallazgos y conclusiones. Formato IEEE.",
        "introduccion": f"Escribe una introduccion academica (400-500 palabras) en espanol sobre {nombre_org}. Incluye: importancia del organismo, contexto biologico, objetivo del estudio (analisis bioinformatico completo del genoma), y estructura del informe. Con referencias a literatura cientifica.",
        "discusion": f"Escribe una discusion academica (400-500 palabras) en espanol sobre los resultados del analisis genomico de {nombre_org}. Discute: significado de las metricas encontradas, comparacion con literatura, implicaciones biologicas, limitaciones del analisis, y perspectivas futuras.",
        "conclusiones": f"Escribe 5-7 conclusiones concisas en espanol para el informe de analisis genomico de {nombre_org}. Cada conclusion en un punto, basada en los datos tipicos de este organismo."
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
    """Limpia texto de IA para uso en ReportLab (quita markdown)."""
    if not texto:
        return ""
    texto = texto.replace("**", "").replace("*", "")
    texto = texto.replace("```", "").replace("`", "")
    texto = texto.replace("##", "").replace("#", "")
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
    elements.append(Paragraph(f"Generado por GenomeHub - {datetime.now().strftime('%d de %B de %Y')}", styles['IEEEAuthor']))
    elements.append(Spacer(1, 0.3 * inch))
    elements.append(Paragraph("Herramientas: BioPython, NCBI GenBank, Gemini AI", styles['IEEEAuthor']))
    elements.append(Paragraph("Formato: IEEE", styles['IEEEAuthor']))
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
            add_body(f"El gen mas largo identificado fue {mayor.get('nombre_gen', 'N/A')} ({mayor.get('producto', 'N/A')}) con {fmt_num(mayor.get('longitud_pb', 0))} pb, mientras que el mas corto fue {menor.get('nombre_gen', 'N/A')} ({menor.get('producto', 'N/A')}) con {fmt_num(menor.get('longitud_pb', 0))} pb.")

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

    # --- 3C. Distancias Intergenicas ---
    dist_key = f"analisis_distancias_{genome}"
    if dist_key in datos:
        dist_data = datos[dist_key]
        add_heading2("Distancias Intergenicas", "C.")

        stats_dist = dist_data.get("estadisticas_generales", {})
        if stats_dist:
            add_body(f"Se analizaron {fmt_num(stats_dist.get('total_pares_adyacentes', 'N/A'))} pares de genes adyacentes. La distancia intergenica promedio fue de {fmt_num(stats_dist.get('distancia_promedio_pb', 0), 1)} pb con una mediana de {fmt_num(stats_dist.get('distancia_mediana_pb', 0), 1)} pb.")
            add_body(f"Se identificaron {fmt_num(stats_dist.get('genes_solapados', 0))} pares de genes solapados ({stats_dist.get('porcentaje_solapados', 0):.1f}% del total) y {stats_dist.get('regiones_grandes_5kb', 0)} regiones intergenicas mayores a 5 kb, que podrian representar islas genomicas o regiones regulatorias.")

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

    # --- 3E. Comparacion de Genomas ---
    for key, data in datos.items():
        if "comparacion_" in key and isinstance(data, dict) and data.get("organismos_comparados"):
            add_heading2("Comparacion de Genomas", "E.")

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

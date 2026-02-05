# Análisis Comparativo del Genoma de E. coli K-12 vs Salmonella LT2

Proyecto de bioinformática para el análisis computacional comparativo entre *Escherichia coli* K-12 MG1655 (comensal) y *Salmonella enterica* serovar Typhimurium LT2 (patógena).

## Descripción

Este proyecto simula el trabajo real de un bioinformático profesional en la industria biotecnológica. Incluye:

- Descarga programática de genomas desde NCBI (E. coli y Salmonella)
- Análisis de codones de inicio (ATG) y parada (TAA, TAG, TGA)
- Extracción y análisis de genes anotados de ambos organismos
- **Análisis comparativo entre comensal y patógeno**
- **Identificación de genes de virulencia**
- **Análisis de distancias intergénicas e islas de patogenicidad**
- 13 visualizaciones profesionales con matplotlib/seaborn
- Validación de resultados con literatura científica
- Documentación reproducible

## Organismos Analizados

| Organismo | Tipo | ID NCBI | Tamaño |
|-----------|------|---------|--------|
| *E. coli* K-12 MG1655 | Comensal (no patógeno) | NC_000913.3 | 4.64 Mb |
| *S. enterica* Typhimurium LT2 | Patógeno | NC_003197.2 | 4.86 Mb |

## Estructura del Proyecto

```
analisis-genoma-ecoli-k12/
│
├── datos/
│   ├── crudo/                         # Datos descargados de NCBI
│   │   ├── ecoli_k12.gb               # E. coli GenBank
│   │   ├── ecoli_k12.fasta            # E. coli FASTA
│   │   ├── salmonella_lt2.gb          # Salmonella GenBank
│   │   ├── salmonella_lt2.fasta       # Salmonella FASTA
│   │   └── metadata_descarga.json
│   └── procesado/                     # Datos procesados
│
├── scripts/
│   ├── descargar_genoma.py            # Descarga desde NCBI (menú interactivo)
│   ├── analisis_codones.py            # Análisis de ATG, TAA, TAG, TGA
│   ├── analisis_genes.py              # Análisis de genes (1=E.coli, 2=Salmonella)
│   ├── comparar_genomas.py            # Comparación E. coli vs Salmonella
│   └── visualizaciones.py             # Generación de 13 gráficos
│
├── resultados/
│   ├── tablas/                        # Archivos JSON y CSV
│   │   ├── analisis_codones_completo.json
│   │   ├── analisis_genes_ecoli_k12.json
│   │   ├── analisis_genes_salmonella_lt2.json
│   │   ├── comparacion_ecoli_vs_salmonella.json
│   │   ├── codones_conteo.csv
│   │   ├── lista_genes_ecoli_k12.csv
│   │   └── lista_genes_salmonella_lt2.csv
│   │
│   └── figuras/                       # 13 Gráficos PNG
│       │
│       │  # Análisis de Codones (3)
│       ├── codones_parada.png
│       ├── codones_inicio_atg.png
│       ├── contenido_gc.png
│       │
│       │  # Análisis de Genes (5)
│       ├── distribucion_tamanos_genes.png
│       ├── distribucion_hebras.png
│       ├── comparacion_literatura.png
│       ├── genes_vs_cds.png
│       ├── genes_extremos.png
│       │
│       │  # Resumen General (1)
│       ├── resumen_general.png
│       │
│       │  # Comparación E. coli vs Salmonella (4)
│       ├── comparacion_genomas_metricas.png
│       ├── comparacion_virulencia.png
│       ├── comparacion_distancias_intergenicas.png
│       └── resumen_comparacion_ecoli_salmonella.png
│
├── informe/
│   └── informe_cientifico.md          # Informe científico final
│
├── requirements.txt                   # Dependencias de Python
└── README.md                          # Este archivo
```

## Requisitos

### Software
- Python 3.8 o superior
- pip (gestor de paquetes de Python)
- Conexión a internet (para descarga de NCBI)

### Librerías de Python
```
biopython>=1.81
pandas>=2.0.0
matplotlib>=3.7.0
seaborn>=0.12.0
numpy>=1.24.0
```

## Instalación

### 1. Clonar el repositorio
```bash
git clone https://github.com/LuisAnthony1/analisis-genoma-ecoli-k12.git
cd analisis-genoma-ecoli-k12
```

### 2. Crear entorno virtual (recomendado)
```bash
python3 -m venv env
source env/bin/activate  # Linux/Mac
# o en Windows: env\Scripts\activate
```

### 3. Instalar dependencias
```bash
pip install -r requirements.txt
```

## Uso

### Pipeline Completo de Análisis

```bash
# 1. Descargar genomas desde NCBI (menú interactivo)
python scripts/descargar_genoma.py
# Seleccionar opción 1 (E. coli) y luego opción 2 (Salmonella)

# 2. Analizar codones de E. coli (ATG, TAA, TAG, TGA)
python scripts/analisis_codones.py

# 3. Analizar genes de E. coli K-12
python scripts/analisis_genes.py 1

# 4. Analizar genes de Salmonella LT2
python scripts/analisis_genes.py 2

# 5. Comparar genomas E. coli vs Salmonella
python scripts/comparar_genomas.py

# 6. Generar todas las visualizaciones (13 gráficos)
python scripts/visualizaciones.py
```

### Ejecución en AWS EC2 (servidor headless)
El script de visualizaciones está optimizado para servidores sin GUI:
- Backend Agg para renderizado sin display
- Gestión de memoria con `gc.collect()`
- DPI optimizado (150) para balance calidad/memoria

## Resultados Principales

### Métricas Comparativas

| Métrica | E. coli K-12 | Salmonella LT2 | Diferencia |
|---------|--------------|----------------|------------|
| Tamaño genoma | 4,641,652 pb | 4,857,450 pb | +215,798 pb (4.65%) |
| Total genes CDS | 4,318 | 4,452 | +134 genes |
| Contenido GC | 50.79% | 52.22% | +1.43% |
| Densidad génica | 87.23% | 86.89% | -0.34% |
| Genes virulencia | 227 | 236 | +9 genes |

### Genes de Virulencia por Categoría

| Categoría | E. coli K-12 | Salmonella LT2 |
|-----------|--------------|----------------|
| Adhesión/Fimbrias | 68 | 51 |
| Sistema de secreción | 21 | 38 |
| Invasión celular | 1 | 22 |
| Captación de hierro | 8 | 22 |
| Resistencia | 11 | 28 |
| Toxinas | 43 | 1 |
| Motilidad/Flagelos | 37 | 33 |

### Hallazgos Clave

1. **Salmonella tiene ~216 kb adicionales** - Este ADN extra contiene genes de virulencia y adaptación al estilo de vida patógeno.

2. **Diferencia en sistemas de secreción** - Salmonella tiene casi el doble de genes de sistemas de secreción tipo III (T3SS), esenciales para inyectar efectores en células huésped.

3. **Invasión celular** - Salmonella tiene 22 genes de invasión vs solo 1 en E. coli, explicando su capacidad de invadir células epiteliales.

4. **Islas de patogenicidad (SPIs)** - Ambos genomas tienen ~8 regiones intergénicas >5kb, pero las de Salmonella contienen SPIs con genes de virulencia.

5. **E. coli K-12 es segura** - Cepa de laboratorio domesticada desde 1922, ha perdido genes de virulencia y es el organismo modelo por excelencia.

## Descripción de Scripts

### `descargar_genoma.py`
- Menú interactivo para seleccionar organismo
- Descarga en formatos GenBank y FASTA
- Usa `gbwithparts` para obtener todas las anotaciones
- Validación de integridad y manejo de errores

### `analisis_codones.py`
- Cuenta codones ATG (inicio) en todo el genoma
- Cuenta codones de parada (TAA, TAG, TGA)
- Calcula densidades por kilobase
- Compara con literatura científica
- Exporta resultados en JSON y CSV

### `analisis_genes.py`
- Argumento: `1` para E. coli, `2` para Salmonella
- Extrae genes CDS del archivo GenBank
- Calcula estadísticas: total, tamaños, GC, distribución por hebras
- Identifica genes más largo y más corto
- Analiza genes codificantes vs no codificantes (tRNA, rRNA, ncRNA)
- Compara con valores de literatura

### `comparar_genomas.py`
- Requiere análisis previo de ambos organismos
- Compara métricas generales de ambos genomas
- Identifica y categoriza genes de virulencia
- Analiza distancias intergénicas e islas genómicas
- Genera resumen interpretativo de diferencias

### `visualizaciones.py`
- Genera 13 gráficos profesionales
- Optimizado para ejecución en servidor (AWS EC2)
- Fallback automático de estilos matplotlib
- Gestión de memoria para entornos limitados

## Visualizaciones Generadas

| # | Archivo | Descripción |
|---|---------|-------------|
| 1 | `codones_parada.png` | Proporciones TAA/TAG/TGA vs literatura |
| 2 | `codones_inicio_atg.png` | Distribución de ATG en el genoma |
| 3 | `contenido_gc.png` | Composición nucleotídica y GC% |
| 4 | `distribucion_tamanos_genes.png` | Histograma de tamaños de genes |
| 5 | `distribucion_hebras.png` | Genes en hebra forward vs reverse |
| 6 | `comparacion_literatura.png` | Valores observados vs literatura |
| 7 | `genes_vs_cds.png` | Genes codificantes vs no codificantes |
| 8 | `genes_extremos.png` | Gen más largo y más corto |
| 9 | `resumen_general.png` | Dashboard de métricas principales |
| 10 | `comparacion_genomas_metricas.png` | E. coli vs Salmonella métricas |
| 11 | `comparacion_virulencia.png` | Genes de virulencia por categoría |
| 12 | `comparacion_distancias_intergenicas.png` | Islas genómicas |
| 13 | `resumen_comparacion_ecoli_salmonella.png` | Resumen comparativo final |

## Validación

Los resultados fueron validados mediante:

1. **Comparación con literatura científica**
   - Blattner et al. (1997) - Secuencia completa de E. coli K-12
   - McClelland et al. (2001) - Secuencia de Salmonella LT2
   - EcoCyc Database - Valores de referencia actualizados

2. **Bases de datos genómicas**
   - NCBI RefSeq NC_000913.3 (E. coli)
   - NCBI RefSeq NC_003197.2 (Salmonella)

3. **Control de calidad**
   - Diferencias <1% con valores de literatura
   - Verificación de rangos esperados

## Referencias

1. Blattner, F. R., et al. (1997). The complete genome sequence of *Escherichia coli* K-12. *Science*, 277(5331), 1453-1462.
2. McClelland, M., et al. (2001). Complete genome sequence of *Salmonella enterica* serovar Typhimurium LT2. *Nature*, 413(6858), 852-856.
3. NCBI Reference Sequences: NC_000913.3, NC_003197.2
4. EcoCyc Database: https://ecocyc.org/

## Autores

| Nombre | Email | Universidad |
|--------|-------|-------------|
| Luis Anthony Mamani Mescco | 194522@unsaac.edu.pe | UNSAAC |
| Tirssa Ivonne Guevara D | 192420@unsaac.edu.pe | UNSAAC |
| Lino Zeynt Huaracallo Arenas | 204798@unsaac.edu.pe | UNSAAC |
| Aracely Llancaya Tapia | 220549@unsaac.edu.pe | UNSAAC |
| Medaly Lozano Llacctahuaman | 195050@unsaac.edu.pe | UNSAAC |

**Universidad**: Universidad Nacional de San Antonio Abad del Cusco (UNSAAC)
**Curso**: Bioinformática
**Fecha**: 2026

## Licencia

Este proyecto es para fines educativos como parte del curso de bioinformática.

---

*Proyecto de análisis genómico comparativo desarrollado como simulación de trabajo profesional en bioinformática.*

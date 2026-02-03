# Análisis del Genoma de E. coli K-12 MG1655

Proyecto de bioinformática para el análisis computacional del genoma completo de la bacteria *Escherichia coli* K-12 MG1655.

## Descripción

Este proyecto simula el trabajo real de un bioinformático profesional en la industria biotecnológica. Incluye:

- Descarga programática del genoma desde NCBI
- Análisis de codones de inicio (ATG) y parada (TAA, TAG, TGA)
- Extracción y análisis de genes anotados
- Visualizaciones profesionales con matplotlib
- Validación de resultados con literatura científica e IA
- Documentación reproducible

## Estructura del Proyecto

```
ecoli_k12_analysis/
│
├── datos/
│   ├── crudo/                    # Datos descargados de NCBI
│   │   ├── ecoli_k12.gb          # Genoma en formato GenBank
│   │   ├── ecoli_k12.fasta       # Genoma en formato FASTA
│   │   └── metadata_descarga.json
│   └── procesado/                # Datos procesados (si aplica)
│
├── scripts/
│   ├── descargar_genoma.py       # Descarga desde NCBI
│   ├── analisis_codones.py       # Análisis de ATG, TAA, TAG, TGA
│   ├── analisis_genes.py         # Extracción y análisis de genes
│   └── visualizaciones.py        # Generación de gráficos
│
├── resultados/
│   ├── tablas/                   # Archivos CSV y JSON
│   │   ├── analisis_codones_completo.json
│   │   ├── analisis_genes_completo.json
│   │   ├── codones_conteo.csv
│   │   ├── lista_genes.csv
│   │   └── estadisticas_genes.csv
│   └── figuras/                  # Gráficos PNG
│       ├── codones_parada.png
│       ├── codones_inicio_atg.png
│       ├── contenido_gc.png
│       ├── distribucion_tamanos_genes.png
│       ├── distribucion_hebras.png
│       ├── comparacion_literatura.png
│       └── resumen_general.png
│
├── documentacion/
│   └── validacion_ia.md          # Validación con IA
│
├── informe/                      # Informe científico final
│
├── requirements.txt              # Dependencias de Python
└── README.md                     # Este archivo
```

## Requisitos

### Software
- Python 3.8 o superior
- pip (gestor de paquetes de Python)

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

### Ejecutar todo el análisis en orden:

```bash
# 1. Descargar el genoma desde NCBI
python scripts/descargar_genoma.py

# 2. Analizar codones (ATG, TAA, TAG, TGA)
python scripts/analisis_codones.py

# 3. Analizar genes anotados
python scripts/analisis_genes.py

# 4. Generar visualizaciones
python scripts/visualizaciones.py
```

## Resultados Esperados

### Métricas del Genoma
| Métrica | Valor Esperado | Fuente |
|---------|----------------|--------|
| Longitud del genoma | ~4,641,652 pb | NCBI |
| Total de genes | ~4,300 | Literatura |
| Contenido GC | ~50.8% | Literatura |
| Densidad génica | ~87% | Literatura |

### Codones de Parada (en genes)
| Codón | Proporción |
|-------|------------|
| TAA | ~63% |
| TGA | ~29% |
| TAG | ~8% |

## Descripción de Scripts

### `descargar_genoma.py`
- Descarga el genoma de E. coli K-12 MG1655 desde NCBI (ID: U00096.3)
- Formatos: GenBank (.gb) y FASTA (.fasta)
- Incluye validación de integridad y manejo de errores

### `analisis_codones.py`
- Cuenta todos los codones ATG en el genoma
- Cuenta codones de parada (TAA, TAG, TGA)
- Calcula densidades por kilobase
- Compara con literatura científica
- Exporta resultados en JSON y CSV

### `analisis_genes.py`
- Extrae genes del archivo GenBank
- Calcula estadísticas: total, tamaños, contenido GC
- Analiza distribución por hebras (forward/reverse)
- Compara con valores de literatura
- Exporta lista completa de genes

### `visualizaciones.py`
- Genera gráficos profesionales con matplotlib/seaborn
- Incluye: histogramas, gráficos circulares, barras comparativas
- Todos los gráficos incluyen títulos y leyendas descriptivas

## Validación

Los resultados fueron validados mediante:
1. **Comparación con literatura científica**: Valores de referencia de Blattner et al. (1997) y bases de datos genómicas
2. **Validación con IA**: Consultas a Claude (Anthropic) para verificar interpretaciones biológicas
3. **Control de calidad**: Verificación de que los valores obtenidos están dentro de rangos esperados

## Referencias

1. Blattner, F. R., et al. (1997). The complete genome sequence of Escherichia coli K-12. Science, 277(5331), 1453-1462.
2. NCBI Reference Sequence: NC_000913.3 / U00096.3
3. EcoGene Database: https://ecogene.org/

## Autor

- **Nombre**: Luis Anthony Mamani Mescco
- **Nombre**: Tirssa Ivonne Guevara D
- **Nombre**: Lino Zeynt Huaracallo Arenas
- **Email**: 194522@unsaac.edu.pe
- **Email**: 192420@unsaac.edu.pe
- **Email**: 204798@unsaac.edu.pe
- **Universidad**: Universidad Nacional de San Antonio Abad del Cusco (UNSAAC)
- **Curso**: BioInformatica
- **Fecha**: 2026

## Licencia

Este proyecto es para fines educativos como parte del curso de bioinformática.

---

*Proyecto desarrollado como simulación de trabajo profesional en bioinformática.*

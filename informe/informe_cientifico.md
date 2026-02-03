# Informe Cientfico: Analisis del Genoma de *Escherichia coli* K-12 MG1655

## Informacion del Proyecto

| Campo | Valor |
|-------|-------|
| **Autor** | Luis Anthony Mamani Mescco|
| | Aracely Llancaya Tapia|
| **Email** | 194522@unsaac.edu.pe |
| **Universidad** | Universidad Nacional de San Antonio Abad del Cusco (UNSAAC) |
| **Curso** | Bioinformatica |
| **Fecha de Analisis** | 2026-02-02 |
| **Servidor** | AWS EC2 Ubuntu (ip-172-31-8-47) |

---

## 1. Introduccion

### 1.1 Contexto Biologico

*Escherichia coli* es una bacteria gram-negativa que habita en el tracto intestinal de mamiferos. La cepa K-12 MG1655 es una de las mas estudiadas en biologia molecular y se considera el organismo modelo por excelencia para estudios geneticos bacterianos.

El genoma de E. coli K-12 fue completamente secuenciado por Blattner et al. en 1997, marcando un hito en la genomica bacteriana. Este genoma sirve como referencia para estudios comparativos y funcionales.

### 1.2 Objetivos

1. Descargar programaticamente el genoma completo desde NCBI
2. Analizar la distribucion de codones de inicio (ATG) y parada (TAA, TAG, TGA)
3. Extraer y caracterizar los genes anotados
4. Comparar los resultados con valores de literatura cientifica
5. Generar visualizaciones profesionales de los datos

---

## 2. Materiales y Metodos

### 2.1 Obtencion de Datos

- **Fuente**: NCBI GenBank
- **Accession Number**: U00096.3
- **Formatos descargados**: GenBank (.gb) y FASTA (.fasta)
- **Herramienta**: BioPython (modulo Entrez)

### 2.2 Herramientas Computacionales

| Herramienta | Version | Proposito |
|-------------|---------|-----------|
| Python | 3.12 | Lenguaje de programacion |
| BioPython | 1.86 | Manejo de secuencias biologicas |
| Pandas | 3.0.0 | Analisis de datos |
| Matplotlib | 3.10.8 | Visualizaciones |
| Seaborn | 0.13.2 | Graficos estadisticos |
| NumPy | 2.4.1 | Calculos numericos |

### 2.3 Infraestructura

- **Servidor**: AWS EC2 t3.micro
- **Sistema Operativo**: Ubuntu 24.04.3 LTS
- **Kernel**: 6.14.0-1018-aws

### 2.4 Scripts Desarrollados

1. `descargar_genoma.py` - Descarga automatizada desde NCBI
2. `analisis_codones.py` - Analisis de ATG, TAA, TAG, TGA
3. `analisis_genes.py` - Extraccion y estadisticas de genes
4. `visualizaciones.py` - Generacion de graficos

---

## 3. Resultados

### 3.1 Caracteristicas Generales del Genoma

| Metrica | Valor Obtenido | Valor Literatura | Diferencia |
|---------|----------------|------------------|------------|
| Longitud del genoma | 4,641,652 pb | 4,641,652 pb | 0 |
| Total de genes (CDS) | 4,318 | ~4,300 | +18 |
| Contenido GC | 50.79% | 50.8% | -0.01% |
| Densidad genica | 87.23% | ~87% | +0.23% |
| Features totales | 9,285 | - | - |

### 3.2 Analisis de Codones de Inicio (ATG)

| Metrica | Valor |
|---------|-------|
| Total de codones ATG | 76,282 |
| Densidad por kilobase | 16.43 ATG/kb |
| Genes anotados | 4,318 |
| Ratio ATG/gen | 17.74 |

**Interpretacion**: Por cada gen anotado existen aproximadamente 18 codones ATG en el genoma. Esto demuestra que no todos los ATG funcionan como codones de inicio reales, ya que:

1. ATG tambien codifica el aminoacido metionina dentro de los genes
2. Muchos ATG aparecen en regiones intergenicas no codificantes
3. Un ATG solo funciona como inicio cuando esta precedido por una secuencia Shine-Dalgarno

### 3.3 Analisis de Codones de Parada

#### 3.3.1 Distribucion en el Genoma Completo

| Codon | Conteo | Densidad/kb | Proporcion |
|-------|--------|-------------|------------|
| TAA | 68,858 | 14.83 | 38.33% |
| TAG | 27,254 | 5.87 | 15.17% |
| TGA | 83,532 | 18.00 | 46.50% |
| **TOTAL** | **179,644** | - | 100% |

#### 3.3.2 Comparacion con Literatura (en genes)

| Codon | Genoma Completo | En Genes (Literatura) | Diferencia |
|-------|-----------------|----------------------|------------|
| TAA | 38.33% | 63% | -24.67% |
| TAG | 15.17% | 8% | +7.17% |
| TGA | 46.50% | 29% | +17.50% |

**Interpretacion**: Las proporciones difieren porque:

- La literatura reporta proporciones **dentro de genes funcionales**
- Nuestro analisis incluye todo el genoma (codificante + no codificante)
- TAA es preferido en genes por mayor eficiencia de terminacion
- Las regiones intergenicas no estan bajo la misma presion selectiva

### 3.4 Contenido de GC

#### 3.4.1 Composicion de Nucleotidos

| Nucleotido | Cantidad | Porcentaje |
|------------|----------|------------|
| Adenina (A) | 1,142,742 | 24.62% |
| Timina (T) | 1,141,382 | 24.59% |
| Guanina (G) | 1,177,437 | 25.37% |
| Citosina (C) | 1,180,091 | 25.42% |

#### 3.4.2 Contenido GC

| Metrica | Valor |
|---------|-------|
| Contenido GC | 50.79% |
| Contenido AT | 49.21% |
| GC en regiones codificantes (CDS) | 50.79% |

**Interpretacion**: El contenido GC de 50.79% indica una composicion equilibrada, tipica de enterobacterias. Este valor:

- Es intermedio comparado con otros organismos
- Afecta la estabilidad del ADN y temperatura de melting
- Influye en el sesgo de uso de codones

### 3.5 Estadisticas de Genes

#### 3.5.1 Tamano de Genes

| Estadistica | Valor |
|-------------|-------|
| Minimo | 27 pb |
| Maximo | 8,622 pb |
| Promedio | 937.64 pb |
| Mediana | 825.00 pb |
| Desviacion Estandar | 649.08 pb |

#### 3.5.2 Distribucion por Tamano

| Categoria | Rango | Cantidad | Porcentaje |
|-----------|-------|----------|------------|
| Muy cortos | 0-300 pb | 532 | 12.32% |
| Cortos | 301-600 pb | 900 | 20.84% |
| Medianos | 601-900 pb | 957 | 22.16% |
| Largos | 901-1500 pb | 1,337 | 30.96% |
| Muy largos | 1501-3000 pb | 533 | 12.34% |
| Extra largos | >3000 pb | 59 | 1.37% |

#### 3.5.3 Distribucion por Hebra

| Hebra | Cantidad | Porcentaje |
|-------|----------|------------|
| Forward (+) | 2,103 | 48.7% |
| Reverse (-) | 2,215 | 51.3% |

**Interpretacion**: La distribucion equilibrada entre hebras (≈50/50) es caracteristica de genomas bacterianos bien organizados y refleja la ausencia de sesgo replicativo significativo.

### 3.6 Tipos de Features en el GenBank

| Tipo | Cantidad |
|------|----------|
| gene | 4,651 |
| CDS | 4,318 |
| ncRNA | 108 |
| tRNA | 86 |
| mobile_element | 50 |
| misc_feature | 48 |
| rRNA | 22 |
| source | 1 |
| rep_origin | 1 |

---

## 4. Visualizaciones Generadas

Se generaron 7 graficos profesionales:

1. **codones_parada.png** - Distribucion de codones de parada
2. **codones_inicio_atg.png** - Analisis de codones ATG
3. **contenido_gc.png** - Composicion de nucleotidos y GC
4. **distribucion_tamanos_genes.png** - Histograma de tamanos de genes
5. **distribucion_hebras.png** - Genes por hebra (forward/reverse)
6. **comparacion_literatura.png** - Valores observados vs literatura
7. **resumen_general.png** - Panel resumen del analisis

---

## 5. Discusion

### 5.1 Validacion de Resultados

Los resultados obtenidos son consistentes con la literatura cientifica:

| Metrica | Concordancia |
|---------|--------------|
| Longitud del genoma | Exacta (4,641,652 pb) |
| Numero de genes | Alta (+18 genes, 0.4% diferencia) |
| Contenido GC | Exacta (50.79% vs 50.8%) |
| Densidad genica | Alta (87.23% vs 87%) |

### 5.2 Diferencias en Codones de Parada

La principal discrepancia observada fue en las proporciones de codones de parada. Esto se explica porque:

1. **Presion selectiva**: Los genes funcionales estan bajo seleccion para usar TAA (mas eficiente)
2. **Regiones no codificantes**: No tienen preferencia por ningun codon de parada
3. **Contexto de medicion**: Literatura mide en genes; nosotros medimos en genoma completo

### 5.3 Significado Biologico

- **Alta densidad genica (87%)**: Tipica de procariotas, refleja genomas compactos y eficientes
- **Distribucion equilibrada por hebra**: Indica ausencia de sesgo replicativo
- **Ratio ATG/gen de 18**: Demuestra la importancia del contexto (Shine-Dalgarno) para la traduccion

---

## 6. Conclusiones

1. Se descargo y analizo exitosamente el genoma completo de E. coli K-12 MG1655
2. Los valores obtenidos coinciden con la literatura cientifica establecida
3. Se identificaron 4,318 genes codificantes con una densidad genica del 87.23%
4. El contenido GC de 50.79% confirma la composicion equilibrada del genoma
5. Las diferencias en proporciones de codones de parada se explican por el contexto de medicion
6. El proyecto demuestra la aplicacion practica de herramientas bioinformaticas en AWS

---

## 7. Archivos Generados

### 7.1 Datos

| Archivo | Ubicacion | Descripcion |
|---------|-----------|-------------|
| ecoli_k12.gb | datos/crudo/ | Genoma en formato GenBank |
| ecoli_k12.fasta | datos/crudo/ | Genoma en formato FASTA |
| metadata_descarga.json | datos/crudo/ | Metadata de la descarga |

### 7.2 Resultados

| Archivo | Ubicacion | Descripcion |
|---------|-----------|-------------|
| analisis_codones_completo.json | resultados/tablas/ | Analisis de codones |
| analisis_genes_completo.json | resultados/tablas/ | Analisis de genes |
| codones_conteo.csv | resultados/tablas/ | Conteo de codones |
| lista_genes.csv | resultados/tablas/ | Lista completa de genes |
| estadisticas_genes.csv | resultados/tablas/ | Estadisticas de genes |

### 7.3 Figuras

Todas las figuras se encuentran en `resultados/figuras/` en formato PNG a 300 DPI.

---

## 8. Referencias

1. Blattner, F. R., et al. (1997). The complete genome sequence of Escherichia coli K-12. *Science*, 277(5331), 1453-1462.

2. NCBI Reference Sequence: NC_000913.3 / U00096.3

3. EcoGene Database: https://ecogene.org/

4. Codon Usage Database (Kazusa): https://www.kazusa.or.jp/codon/

5. BioPython Documentation: https://biopython.org/

---

## 9. Anexo: Comandos Ejecutados

```bash
# Activar entorno virtual
source ~/bioinfo_env/bin/activate

# Instalar dependencias
pip install -r ~/ecoli_k12_analysis/requirements.txt

# Ejecutar scripts
cd ~/ecoli_k12_analysis/scripts
python descargar_genoma.py
python analisis_codones.py
python analisis_genes.py
python visualizaciones.py
```

---

*Informe generado como parte del proyecto de Bioinformatica - UNSAAC 2026*

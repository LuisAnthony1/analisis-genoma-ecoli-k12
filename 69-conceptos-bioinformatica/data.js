// data.js
// Datos para el Mapa Mental 3D de Bioinformática
// Organizados en 9 bloques biológicos con colores del README

// Bloques biológicos
const BLOQUES = [
  {
    id: "bloque1",
    nombre: "Fundamentos Celulares",
    color: "#10b981" // emerald
  },
  {
    id: "bloque2",
    nombre: "Macromoléculas Informacionales",
    color: "#06b6d4" // cyan
  },
  {
    id: "bloque3",
    nombre: "Organización Genómica",
    color: "#a855f7" // purple
  },
  {
    id: "bloque4",
    nombre: "Células y Organización",
    color: "#f59e0b" // amber
  },
  {
    id: "bloque5",
    nombre: "Flujo de Información Genética",
    color: "#ef4444" // red
  },
  {
    id: "bloque6",
    nombre: "Tecnologías ADN Recombinante",
    color: "#ec4899" // pink
  },
  {
    id: "bloque7",
    nombre: "Técnicas de Análisis",
    color: "#3b82f6" // blue
  },
  {
    id: "bloque8",
    nombre: "Variación y Evolución",
    color: "#64748b" // slate
  },
  {
    id: "bloque9",
    nombre: "Bioinformática",
    color: "#f97316" // orange
  }
];

// Nodos: 69 conceptos con bloque y descripcion enriquecida
const NODOS = [
  // === Bloque 1: Fundamentos Celulares (4) ===
  {
    id: "celula",
    nombre: "Célula",
    bloque: "bloque1",
    descripcion:
      "Unidad estructural y funcional basica de todos los organismos vivos. Delimitada por una membrana plasmatica que regula el intercambio de sustancias con el medio externo. Contiene el material genetico (ADN) y la maquinaria molecular necesaria para el metabolismo, la replicacion y la respuesta al entorno.\n\nExisten dos tipos fundamentales: procariotas (sin nucleo definido, como bacterias) y eucariotas (con nucleo, como celulas animales y vegetales). El tamano tipico varia de 1-10 um en procariotas a 10-100 um en eucariotas.\n\nEn bioinformatica, la celula es el contexto biologico donde operan todos los genes, proteinas y vias metabolicas que se analizan computacionalmente. Los modelos de celula completa (whole-cell models) intentan simular todos sus procesos."
  },
  {
    id: "teoria_celular",
    nombre: "Teoría celular",
    bloque: "bloque1",
    descripcion:
      "Principio fundamental de la biologia formulado por Schleiden, Schwann y Virchow en el siglo XIX. Establece tres postulados:\n\n1. Todos los seres vivos estan compuestos por una o mas celulas.\n2. La celula es la unidad basica de estructura y funcion.\n3. Toda celula proviene de una celula preexistente (omnis cellula e cellula).\n\nEste principio implica que la informacion genetica se transmite de celula madre a celulas hijas, lo cual es la base de la herencia. La teoria celular es el cimiento sobre el que se construyen la genetica, la genomica y, por extension, la bioinformatica moderna."
  },
  {
    id: "ciclo_vida_celular",
    nombre: "Ciclo de vida celular",
    bloque: "bloque1",
    descripcion:
      "Serie ordenada de eventos que una celula atraviesa desde su formacion hasta su division en dos celulas hijas. Se compone de:\n\n- Fase G1: crecimiento celular y sintesis de proteinas.\n- Fase S: replicacion del ADN (duplicacion del genoma completo).\n- Fase G2: preparacion para la division.\n- Fase M (mitosis): division del nucleo y citoplasma.\n\nEn bacterias como E. coli, el ciclo puede completarse en tan solo 20 minutos bajo condiciones optimas. Puntos de control (checkpoints) verifican la integridad del ADN en cada transicion.\n\nEn bioinformatica, el analisis de expresion genica durante el ciclo celular revela que genes se activan en cada fase, permitiendo entender la regulacion del crecimiento y la division."
  },
  {
    id: "via_metabolica",
    nombre: "Vía metabólica / Pathway",
    bloque: "bloque1",
    descripcion:
      "Serie de reacciones quimicas conectadas dentro de una celula, donde el producto de una reaccion es el sustrato de la siguiente. Cada paso esta catalizado por una enzima especifica codificada por un gen.\n\nSe clasifican en:\n- Catabolicas: degradan moleculas para obtener energia (ej: glucolisis, ciclo de Krebs).\n- Anabolicas: sintetizan moleculas complejas a partir de precursores simples (ej: sintesis de aminoacidos).\n- De senalizacion: transmiten senales intra/extracelulares.\n\nBases de datos como KEGG y MetaCyc contienen miles de vias metabolicas mapeadas. En bioinformatica, la reconstruccion metabolica a partir del genoma permite predecir las capacidades metabolicas de un organismo, lo cual es esencial en biotecnologia y medicina personalizada."
  },

  // === Bloque 2: Macromoleculas Informacionales (9) ===
  {
    id: "adn",
    nombre: "ADN",
    bloque: "bloque2",
    descripcion:
      "Acido desoxirribonucleico: polimero de nucleotidos que almacena la informacion genetica hereditaria. Compuesto por cuatro bases: adenina (A), timina (T), guanina (G) y citosina (C), unidas a un esqueleto de azucar desoxirribosa y fosfato.\n\nCaracteristicas clave:\n- Doble hebra antiparalela (5'→3' y 3'→5').\n- Estable quimicamente gracias al azucar desoxirribosa.\n- En E. coli K-12: cromosoma circular de ~4.6 millones de pares de bases.\n- En humanos: ~3,200 millones de pb distribuidos en 23 pares de cromosomas.\n\nEl ADN es la molecula central de la bioinformatica. Su secuenciacion, ensamblaje, anotacion y comparacion son las operaciones fundamentales del analisis genomico. Los formatos FASTA y GenBank son los estandares para representar secuencias de ADN digitalmente."
  },
  {
    id: "arn",
    nombre: "ARN",
    bloque: "bloque2",
    descripcion:
      "Acido ribonucleico: polimero de nucleotidos que participa en la transferencia y procesamiento de informacion genetica. Usa uracilo (U) en lugar de timina (T) y ribosa como azucar.\n\nTipos principales:\n- ARNm (mensajero): transporta la informacion del gen al ribosoma.\n- ARNt (transferencia): adapta codones a aminoacidos durante la traduccion.\n- ARNr (ribosomal): componente estructural y catalitico del ribosoma.\n- ARN no codificante: miARN, siARN, lncARN con funciones regulatorias.\n\nGeneralmente monocatenario, puede formar estructuras secundarias complejas (horquillas, bucles). En bioinformatica, el analisis de RNA-seq permite cuantificar la expresion de todos los genes de un organismo simultaneamente, revelando que genes estan activos y en que condiciones."
  },
  {
    id: "proteina",
    nombre: "Proteína",
    bloque: "bloque2",
    descripcion:
      "Macromolecula formada por una o mas cadenas de aminoacidos (polipeptidos) que se pliegan en estructuras tridimensionales especificas para ejecutar funciones biologicas.\n\nFunciones principales:\n- Enzimas: catalizan reacciones quimicas (ej: ADN polimerasa).\n- Estructurales: dan forma a celulas y tejidos (ej: collageno).\n- Transporte: mueven moleculas (ej: hemoglobina transporta O2).\n- Senalizacion: hormonas y receptores.\n- Defensa: anticuerpos del sistema inmune.\n\nLa secuencia de aminoacidos determina la estructura 3D y la funcion. En E. coli K-12, se han identificado ~4,300 genes codificantes de proteinas. Herramientas como AlphaFold predicen la estructura 3D a partir de la secuencia, y bases de datos como UniProt almacenan millones de secuencias anotadas."
  },
  {
    id: "nucleotido",
    nombre: "Nucleótido",
    bloque: "bloque2",
    descripcion:
      "Unidad monomera de los acidos nucleicos (ADN y ARN). Cada nucleotido esta compuesto por tres partes:\n\n1. Base nitrogenada: purina (A, G) o pirimidina (C, T/U) - porta la informacion.\n2. Azucar pentosa: desoxirribosa (ADN) o ribosa (ARN).\n3. Grupo fosfato: confiere carga negativa y forma el enlace fosfodiester entre nucleotidos.\n\nLos nucleotidos se unen por enlaces fosfodiester 5'→3' formando cadenas. Tambien funcionan como moleculas de senalizacion (cAMP, cGMP) y como moneda energetica (ATP = adenosin trifosfato, la principal fuente de energia celular).\n\nEn bioinformatica, cada nucleotido se representa por una letra (A, T/U, G, C), formando las secuencias que se almacenan en bases de datos como GenBank."
  },
  {
    id: "bases_nitrogenadas",
    nombre: "Bases nitrogenadas",
    bloque: "bloque2",
    descripcion:
      "Moleculas aromaticas que portan la informacion genetica en los acidos nucleicos. Se clasifican en:\n\nPurinas (dos anillos fusionados):\n- Adenina (A): se aparea con T (ADN) o U (ARN).\n- Guanina (G): se aparea con C.\n\nPirimidinas (un solo anillo):\n- Citosina (C): se aparea con G (3 puentes de hidrogeno).\n- Timina (T): exclusiva del ADN, se aparea con A (2 puentes de hidrogeno).\n- Uracilo (U): reemplaza a T en el ARN.\n\nEl contenido de GC (porcentaje de G+C) varia entre organismos: E. coli tiene ~50.8% GC, mientras que Plasmodium falciparum tiene solo ~19%. El contenido GC afecta la estabilidad termica del ADN y es un parametro fundamental en bioinformatica para la caracterizacion genomica y el diseno de cebadores para PCR."
  },
  {
    id: "regla_chargaff",
    nombre: "Regla de Chargaff",
    bloque: "bloque2",
    descripcion:
      "Observacion empirica descubierta por Erwin Chargaff en 1950 que establece que en cualquier molecula de ADN bicatenario:\n\n- La cantidad de adenina es aproximadamente igual a la de timina (A ≈ T).\n- La cantidad de guanina es aproximadamente igual a la de citosina (G ≈ C).\n- Por tanto: purinas (A+G) ≈ pirimidinas (C+T).\n\nEsta regla fue clave para que Watson y Crick propusieran el modelo de la doble helice en 1953, ya que implica un apareamiento especifico entre bases complementarias.\n\nEn bioinformatica, la regla de Chargaff se usa como control de calidad: si un genoma secuenciado no cumple estas proporciones en su hebra doble, puede indicar errores de secuenciacion o ensamblaje."
  },
  {
    id: "doble_helice",
    nombre: "Doble hélice",
    bloque: "bloque2",
    descripcion:
      "Estructura tridimensional del ADN descubierta por James Watson y Francis Crick en 1953 (con datos de Rosalind Franklin). Consiste en dos hebras de nucleotidos enrolladas una alrededor de la otra en sentido antiparalelo.\n\nCaracteristicas:\n- Hebras unidas por puentes de hidrogeno entre bases complementarias (A-T: 2 puentes, G-C: 3 puentes).\n- El surco mayor y el surco menor permiten la interaccion con proteinas reguladoras.\n- Giro completo cada ~10.5 pares de bases (forma B, la mas comun).\n- Diametro constante de ~2 nm.\n\nLa forma B es la mas frecuente in vivo, pero existen formas A (deshidratada) y Z (giro a la izquierda). La estructura de doble helice es la que permite la replicacion semiconservativa y la reparacion del ADN, ya que cada hebra contiene la informacion para reconstruir la otra."
  },
  {
    id: "complementariedad_bases",
    nombre: "Complementariedad de bases",
    bloque: "bloque2",
    descripcion:
      "Principio fundamental por el cual las bases nitrogenadas se aparean de forma especifica:\n\n- A se aparea con T (en ADN) o U (en ARN) mediante 2 puentes de hidrogeno.\n- G se aparea con C mediante 3 puentes de hidrogeno.\n\nEste principio es la base de:\n- Replicacion del ADN: cada hebra sirve de molde para sintetizar la complementaria.\n- Transcripcion: la ARN polimerasa lee una hebra y sintetiza ARN complementario.\n- Hibridacion: union de secuencias complementarias de distintas moleculas.\n- PCR: los cebadores se unen por complementariedad a sus sitios diana.\n- Microarreglos: sondas hibridan con transcritos complementarios.\n\nEn bioinformatica, la complementariedad permite calcular la hebra reversa-complementaria de cualquier secuencia, predecir estructuras secundarias de ARN y disenar oligonucleotidos (primers, sondas)."
  },
  {
    id: "replicacion_adn",
    nombre: "Replicación del ADN",
    bloque: "bloque2",
    descripcion:
      "Proceso mediante el cual la celula copia fielmente su ADN antes de dividirse. Es semiconservativa: cada molecula hija contiene una hebra original y una nueva.\n\nPasos principales:\n1. Helicasa desenrolla la doble helice separando las hebras.\n2. Primasa sintetiza un cebador de ARN.\n3. ADN polimerasa III sintetiza la nueva hebra en direccion 5'→3'.\n4. Hebra lider se sintetiza continuamente; hebra rezagada en fragmentos de Okazaki.\n5. ADN ligasa une los fragmentos.\n\nEn E. coli, la replicacion comienza en un unico origen (oriC) y avanza bidireccionalmente a ~1,000 nucleotidos/segundo. La tasa de error es de ~1 en 10^9 bases gracias a la actividad correctora (proofreading) de la polimerasa.\n\nErrores no corregidos durante la replicacion generan mutaciones, que son la materia prima de la evolucion."
  },

  // === Bloque 3: Organizacion Genomica (7) ===
  {
    id: "cromosoma",
    nombre: "Cromosoma",
    bloque: "bloque3",
    descripcion:
      "Estructura organizada de ADN y proteinas que contiene parte o la totalidad del material genetico de un organismo.\n\nEn procariotas:\n- Generalmente un cromosoma circular unico (ej: E. coli tiene 1 cromosoma de ~4.6 Mb).\n- Puede tener plasmidos adicionales (ADN circular extracromosomico).\n- ADN asociado a proteinas tipo histona (nucleoid-associated proteins).\n\nEn eucariotas:\n- Multiples cromosomas lineales (humanos: 46 cromosomas).\n- ADN empaquetado con histonas formando cromatina.\n- Regiones especializadas: centromero (division), telomeros (extremos protectores).\n\nEn bioinformatica, cada cromosoma se representa como una secuencia continua (o scaffolds/contigs si el ensamblaje esta incompleto). El cariotipo y la organizacion cromosomica se analizan mediante genomica comparativa."
  },
  {
    id: "gen",
    nombre: "Gen",
    bloque: "bloque3",
    descripcion:
      "Segmento de ADN que contiene la informacion necesaria para producir un producto funcional (proteina o ARN). Un gen incluye:\n\n- Region promotora: donde se une la ARN polimerasa para iniciar la transcripcion.\n- Region codificante: secuencia que se transcribe (y en genes proteicos, se traduce).\n- Region terminadora: senal de finalizacion de la transcripcion.\n- En eucariotas: exones (codificantes) e intrones (no codificantes, removidos por splicing).\n\nE. coli K-12 tiene ~4,300 genes codificantes de proteinas en ~4.6 Mb. En humanos, solo ~20,000-25,000 genes codifican proteinas en ~3,200 Mb (el resto es ADN intergénico, repetitivo, regulatorio, etc.).\n\nLa anotacion genica (gene annotation) es una tarea central de la bioinformatica: identificar donde estan los genes, que codifican y como se regulan."
  },
  {
    id: "herencia_mendeliana",
    nombre: "Herencia mendeliana",
    bloque: "bloque3",
    descripcion:
      "Patrones de transmision de rasgos descubiertos por Gregor Mendel (1866) trabajando con guisantes. Establece:\n\nLey de segregacion: cada individuo tiene dos alelos por gen, y cada gameto recibe solo uno.\nLey de distribucion independiente: genes en cromosomas distintos se heredan independientemente.\nDominancia: un alelo puede enmascarar la expresion de otro.\n\nTerminologia clave:\n- Alelo: variante de un gen.\n- Homocigoto: dos alelos iguales (AA o aa).\n- Heterocigoto: dos alelos diferentes (Aa).\n- Fenotipo: rasgo observable.\n- Genotipo: composicion genetica.\n\nAunque Mendel no conocia el ADN, sus leyes se explican por la meiosis y la estructura cromosomica. En bioinformatica, los analisis de GWAS (Genome-Wide Association Studies) buscan variantes geneticas asociadas a fenotipos usando principios mendelianos a escala genomica."
  },
  {
    id: "mutacion",
    nombre: "Mutación",
    bloque: "bloque3",
    descripcion:
      "Cambio permanente en la secuencia de nucleotidos del ADN. Las mutaciones son la fuente primaria de variacion genetica.\n\nTipos principales:\n- Puntuales: sustitucion de una base (transicion: purina↔purina; transversion: purina↔pirimidina).\n- Inserciones/deleciones (indels): adicion o eliminacion de nucleotidos.\n- Frameshift: indels que no son multiplo de 3, alterando el marco de lectura.\n\nEfectos en proteinas:\n- Sinonima (silenciosa): no cambia el aminoacido (por degeneracion del codigo).\n- No sinonima (missense): cambia un aminoacido.\n- Sin sentido (nonsense): genera un codon de parada prematuro.\n\nTasa de mutacion en E. coli: ~5.4 × 10^-10 por par de bases por generacion. Las mutaciones pueden ser neutras, daninas o beneficiosas. En bioinformatica, la deteccion de variantes (variant calling) compara secuencias contra un genoma de referencia para identificar mutaciones."
  },
  {
    id: "ligamiento_genetico",
    nombre: "Ligamiento genético",
    bloque: "bloque3",
    descripcion:
      "Fenomeno por el cual genes localizados cerca uno del otro en el mismo cromosoma tienden a heredarse juntos, ya que es menos probable que un evento de recombinacion (crossing-over) los separe.\n\nConceptos clave:\n- Genes completamente ligados: siempre se heredan juntos (recombinacion = 0%).\n- Genes parcialmente ligados: se recombinan con frecuencia < 50%.\n- Genes independientes: en cromosomas distintos, recombinacion = 50%.\n\nLa frecuencia de recombinacion se mide en centiMorgans (cM): 1 cM ≈ 1% de probabilidad de recombinacion por meiosis. En E. coli, la transferencia de genes por conjugacion permitio construir los primeros mapas geneticos bacterianos.\n\nEn bioinformatica, el analisis de ligamiento se usa para mapear genes asociados a enfermedades y para detectar desequilibrio de ligamiento (LD) en estudios de asociacion genomica."
  },
  {
    id: "mapa_genetico",
    nombre: "Mapa genético",
    bloque: "bloque3",
    descripcion:
      "Representacion del orden relativo de genes o marcadores en un cromosoma, basada en las frecuencias de recombinacion entre ellos.\n\nTipos:\n- Mapa de ligamiento: posiciones relativas basadas en frecuencia de recombinacion (unidades: centiMorgans).\n- Mapa fisico: posiciones absolutas en pares de bases sobre el ADN.\n- Mapa citogenetico: posiciones visibles al microscopio en cromosomas tenidos.\n\nHistoria: el primer mapa genetico fue construido por Alfred Sturtevant en 1913 para Drosophila melanogaster. El mapa genetico de E. coli fue pionero, mapeando genes por transferencia interrumpida en conjugacion (experimentos de Jacob y Wollman).\n\nEn bioinformatica moderna, los mapas geneticos se generan computacionalmente alineando marcadores moleculares (SNPs, microsatelites) a genomas de referencia."
  },
  {
    id: "un_gen_una_proteina",
    nombre: "Un gen - una proteina",
    bloque: "bloque3",
    descripcion:
      "Hipotesis formulada originalmente como 'un gen–una enzima' por Beadle y Tatum en 1941, basandose en mutantes de Neurospora crassa. Propone que cada gen codifica una sola enzima.\n\nFue refinada a 'un gen–una proteina' al descubrirse proteinas no enzimaticas, y luego a 'un gen–un polipeptido' al reconocer que proteinas pueden tener multiples subunidades.\n\nHoy sabemos que esta hipotesis es una simplificacion:\n- El splicing alternativo permite que un gen produzca multiples proteinas.\n- Un gen puede codificar ARN no codificante funcional.\n- Las modificaciones post-traduccionales diversifican aun mas los productos.\n\nEn humanos, ~20,000 genes producen >100,000 proteinas distintas gracias al splicing alternativo. Sin embargo, en bacterias como E. coli, la correspondencia un gen–una proteina se mantiene bastante bien, ya que carecen de intrones."
  },

  // === Bloque 4: Celulas y Organizacion (4) ===
  {
    id: "celula_eucariota",
    nombre: "Célula eucariota",
    bloque: "bloque4",
    descripcion:
      "Celula que posee un nucleo definido rodeado por membrana nuclear, donde se almacena el ADN organizado en cromosomas lineales.\n\nCaracteristicas distintivas:\n- Organulos membranosos: nucleo, mitocondrias, reticulo endoplasmatico, aparato de Golgi.\n- Genoma mas grande y complejo que procariotas (Mb a Gb).\n- Genes interrumpidos por intrones que requieren splicing.\n- Transcripcion en el nucleo, traduccion en el citoplasma (separacion espacial y temporal).\n- Regulacion compleja: epigenetica (metilacion, modificacion de histonas), enhancers distales.\n\nEjemplos: levaduras (Saccharomyces cerevisiae, ~12 Mb), celulas humanas (~3,200 Mb). En bioinformatica, el analisis de genomas eucariotas requiere herramientas especializadas para predecir intrones/exones, promotores y elementos regulatorios distales."
  },
  {
    id: "celula_procariota",
    nombre: "Célula procariota",
    bloque: "bloque4",
    descripcion:
      "Celula sin nucleo definido ni organulos membranosos. El ADN se encuentra en una region del citoplasma llamada nucleoide.\n\nCaracteristicas distintivas:\n- Cromosoma principal generalmente circular, mas plasmidos opcionales.\n- Genes continuos (sin intrones, salvo raras excepciones).\n- Operones: grupos de genes co-transcritos bajo un unico promotor.\n- Transcripcion y traduccion acopladas: el ribosoma traduce el ARNm mientras aun se esta transcribiendo.\n- Genomas compactos: E. coli K-12 tiene ~4.6 Mb con ~87% de secuencia codificante.\n\nIncluyen bacterias y arqueas. Las bacterias son el principal objeto de estudio en bioinformatica genomica por sus genomas pequenos, bien anotados y abundantes. El Proyecto Genoma de E. coli fue completado en 1997, siendo uno de los primeros genomas bacterianos secuenciados."
  },
  {
    id: "exon",
    nombre: "Exón",
    bloque: "bloque4",
    descripcion:
      "Segmento de un gen que se retiene en el ARN mensajero maduro despues del splicing y que tipicamente contribuye a la secuencia codificante de la proteina.\n\nCaracteristicas:\n- En eucariotas, un gen tipico tiene multiples exones separados por intrones.\n- El primer exon contiene la secuencia 5'UTR y el codon de inicio.\n- El ultimo exon contiene el codon de parada y la senal de poliadenilacion.\n- El splicing alternativo puede incluir o excluir ciertos exones, generando isoformas proteicas diferentes.\n\nEjemplo: el gen de la titina humana tiene 363 exones, produciendo la proteina mas grande conocida (34,350 aminoacidos).\n\nEn bioinformatica, la prediccion de exones es una tarea fundamental de la anotacion genomica. Programas como Augustus, GeneMark y BRAKER utilizan modelos de Markov ocultos (HMM) y datos de RNA-seq para identificar las fronteras exon-intron."
  },
  {
    id: "intron",
    nombre: "Intrón",
    bloque: "bloque4",
    descripcion:
      "Segmento interveniente dentro de un gen que es transcrito en ARN pero eliminado (cortado) durante el splicing, no apareciendo en el ARNm maduro.\n\nCaracteristicas:\n- Presentes principalmente en eucariotas; muy raros en procariotas.\n- Contienen secuencias conservadas en sus bordes: GU en el extremo 5' (donor) y AG en el extremo 3' (acceptor).\n- El espliceosoma (splicing complex) es la maquinaria que los remueve.\n- Pueden contener elementos regulatorios (enhancers intronicos).\n\nDatos curiosos:\n- En humanos, los intrones representan ~25% del genoma total.\n- Algunos genes tienen intrones de >100 kb (mayores que genomas bacterianos completos).\n- En E. coli, los genes carecen de intrones, lo que simplifica la anotacion genomica.\n\nEn bioinformatica, la deteccion de intrones es critica para la anotacion correcta de genes eucariotas. RNA-seq mapea lecturas que cruzan fronteras exon-exon (junction reads), confirmando la presencia de intrones."
  },

  // === Bloque 5: Flujo de Informacion Genetica (16) ===
  {
    id: "transcripcion",
    nombre: "Transcripción",
    bloque: "bloque5",
    descripcion:
      "Primer paso del flujo de informacion genetica: sintesis de una molecula de ARN usando una hebra de ADN como molde. La ARN polimerasa lee la hebra molde en direccion 3'→5' y sintetiza el ARN en direccion 5'→3'.\n\nEtapas:\n1. Iniciacion: la ARN polimerasa reconoce y se une al promotor del gen.\n2. Elongacion: la enzima avanza sintetizando el ARN complementario.\n3. Terminacion: la enzima se detiene al llegar a una senal de terminacion.\n\nEn procariotas: una sola ARN polimerasa transcribe todos los tipos de genes. Los operones permiten transcribir multiples genes en un solo ARNm (policistronico).\n\nEn eucariotas: tres ARN polimerasas (I, II, III) para distintos tipos de ARN. ARN pol II transcribe genes codificantes de proteinas.\n\nEn bioinformatica, RNA-seq mide los niveles de transcripcion de todos los genes simultaneamente, permitiendo estudiar la expresion diferencial entre condiciones."
  },
  {
    id: "arn_polimerasa",
    nombre: "ARN polimerasa",
    bloque: "bloque5",
    descripcion:
      "Complejo enzimatico multi-subunidad que cataliza la sintesis de ARN a partir de un molde de ADN, incorporando ribonucleotidos complementarios.\n\nEn E. coli:\n- Core enzyme: subunidades α2ββ'ω.\n- Holoenzima: core + factor sigma (σ70 para genes housekeeping).\n- Velocidad: ~40-80 nucleotidos/segundo.\n- No requiere cebador (a diferencia de la ADN polimerasa).\n\nEn eucariotas existen tres tipos:\n- ARN Pol I: sintetiza ARNr (45S) en el nucleolo.\n- ARN Pol II: sintetiza ARNm y algunos ARN pequenos.\n- ARN Pol III: sintetiza ARNt, ARNr 5S y otros ARN pequenos.\n\nLa ARN polimerasa es diana de antibioticos como la rifampicina (inhibe la ARN polimerasa bacteriana, usada contra tuberculosis). En bioinformatica, los datos de ChIP-seq de ARN polimerasa revelan que genes estan siendo activamente transcritos."
  },
  {
    id: "arn_mensajero",
    nombre: "ARN mensajero (ARNm)",
    bloque: "bloque5",
    descripcion:
      "Molecula de ARN que transporta la informacion codificante de un gen desde el ADN hasta el ribosoma, donde se traduce en proteina.\n\nEstructura del ARNm eucariota maduro:\n- Cap 5' (7-metilguanosina): protege de degradacion e inicia traduccion.\n- 5'UTR: region no traducida que contiene elementos regulatorios.\n- ORF (Marco de Lectura Abierto): secuencia codificante desde AUG hasta codon stop.\n- 3'UTR: region no traducida con senales de estabilidad y localizacion.\n- Cola poli-A: ~200 adeninas que protegen contra degradacion.\n\nEn procariotas:\n- Sin cap ni cola poli-A.\n- ARNm policistronico: un transcrito codifica multiples proteinas.\n- Secuencia Shine-Dalgarno recluta al ribosoma.\n\nVida media: segundos a minutos en bacterias; horas a dias en eucariotas. En bioinformatica, la prediccion de ORFs es fundamental para la anotacion de genomas."
  },
  {
    id: "splicing",
    nombre: "Splicing (empalme)",
    bloque: "bloque5",
    descripcion:
      "Proceso de maduracion del ARN pre-mensajero eucariota en el que se eliminan los intrones y se unen los exones para formar el ARNm maduro.\n\nMecanismo:\n1. El espliceosoma (complejo de 5 snRNP: U1, U2, U4, U5, U6) reconoce las secuencias consenso en los bordes del intron.\n2. Se forma un lazo (lariat) mediante un ataque nucleofilico de una A del branch point.\n3. El intron se libera como lariat y los exones se unen.\n\nSplicing alternativo: permite que un mismo gen produzca diferentes ARNm al incluir/excluir exones selectivamente. En humanos, ~95% de los genes multiexonicos sufren splicing alternativo.\n\nEn procariotas: no hay splicing convencional (no hay intrones en la mayoria de genes), lo que simplifica enormemente la prediccion de genes. En bioinformatica, herramientas como STAR y HISAT2 alinean lecturas de RNA-seq considerando junctions de splicing."
  },
  {
    id: "traduccion",
    nombre: "Traducción",
    bloque: "bloque5",
    descripcion:
      "Segundo paso principal del flujo de informacion genetica: sintesis de una proteina usando la informacion codificada en el ARNm. El ribosoma lee codones (tripletes) y los ARNt aportan los aminoacidos correspondientes.\n\nEtapas:\n1. Iniciacion: la subunidad menor del ribosoma se une al ARNm, reconoce el AUG de inicio, y se une la subunidad mayor.\n2. Elongacion: el ribosoma avanza codon por codon, incorporando aminoacidos mediante enlaces peptidicos.\n3. Terminacion: un codon stop (UAA, UAG, UGA) detiene la traduccion; factores de liberacion separan el polipeptido.\n\nVelocidad: ~20 aminoacidos/segundo en E. coli.\nMultiples ribosomas pueden traducir simultaneamente un mismo ARNm (polisomas).\n\nEn bioinformatica, el analisis de uso de codones revela preferencias de cada organismo (codon usage bias), y la proteómica complementa la genómica al medir los productos finales de la traduccion."
  },
  {
    id: "ribosoma",
    nombre: "Ribosoma",
    bloque: "bloque5",
    descripcion:
      "Complejo macromolecular formado por ARN ribosomal (ARNr) y proteinas que cataliza la sintesis de proteinas (traduccion).\n\nEstructura:\n- Procariotas (70S): subunidad menor 30S (ARNr 16S + 21 proteinas) + subunidad mayor 50S (ARNr 23S, 5S + 31 proteinas).\n- Eucariotas (80S): subunidad menor 40S (ARNr 18S) + subunidad mayor 60S (ARNr 28S, 5.8S, 5S).\n\nSitios funcionales:\n- Sitio A (aminoacil): donde entra el ARNt cargado.\n- Sitio P (peptidil): donde esta el ARNt con la cadena en crecimiento.\n- Sitio E (exit): por donde sale el ARNt descargado.\n\nEl ARNr 16S es la molecula mas usada en filogenetica bacteriana: su secuencia se usa para clasificar y construir arboles evolutivos. Es la base del proyecto de microbiomas como el Human Microbiome Project. Antibioticos como tetraciclina, eritromicina y cloranfenicol actuan inhibiendo el ribosoma bacteriano."
  },
  {
    id: "aminoacido",
    nombre: "Aminoácido",
    bloque: "bloque5",
    descripcion:
      "Monomero de las proteinas. Cada aminoacido tiene una estructura comun: un grupo amino (-NH2), un grupo carboxilo (-COOH), un hidrogeno y una cadena lateral (R) que determina sus propiedades.\n\n20 aminoacidos estandar codificados por el codigo genetico:\n- No polares: Gly, Ala, Val, Leu, Ile, Pro, Phe, Met, Trp.\n- Polares sin carga: Ser, Thr, Cys, Tyr, Asn, Gln.\n- Cargados positivamente: Lys, Arg, His.\n- Cargados negativamente: Asp, Glu.\n\nPropiedades importantes:\n- La cadena lateral determina la hidrofobicidad, carga y funcion.\n- La cisteina puede formar puentes disulfuro (estabilidad estructural).\n- Aminoacidos esenciales deben obtenerse de la dieta en humanos.\n\nEn bioinformatica, las matrices de sustitucion (BLOSUM, PAM) cuantifican la probabilidad de reemplazo entre aminoacidos durante la evolucion, y son la base de los alineamientos de secuencias proteicas."
  },
  {
    id: "polipeptido",
    nombre: "Polipéptido",
    bloque: "bloque5",
    descripcion:
      "Cadena lineal de aminoacidos unidos por enlaces peptidicos (enlace covalente entre el grupo carboxilo de un aminoacido y el grupo amino del siguiente, con liberacion de agua).\n\nNiveles de organizacion:\n- Estructura primaria: secuencia lineal de aminoacidos.\n- Estructura secundaria: patrones locales (alfa-helice, lamina-beta, giros).\n- Estructura terciaria: plegamiento 3D completo de la cadena.\n- Estructura cuaternaria: ensamblaje de multiples cadenas polipeptidicas.\n\nUn polipeptido se convierte en proteina funcional al plegarse correctamente y, en muchos casos, al recibir modificaciones post-traduccionales (fosforilacion, glicosilacion, acetilacion, etc.).\n\nEn bioinformatica, herramientas como AlphaFold y ESMFold predicen la estructura 3D a partir de la secuencia de aminoacidos con precision atomica, revolucionando la biologia estructural."
  },
  {
    id: "codon",
    nombre: "Codón",
    bloque: "bloque5",
    descripcion:
      "Secuencia de tres nucleotidos consecutivos en el ARNm que especifica un aminoacido o una senal de control.\n\nCaracteristicas:\n- Hay 64 codones posibles (4^3 = 64 combinaciones de bases).\n- 61 codifican aminoacidos; 3 son codones de parada (UAA, UAG, UGA).\n- AUG codifica metionina Y sirve como codon de inicio.\n- La lectura es sin solapamiento y sin espacios (marco de lectura continuo).\n\nEl marco de lectura (reading frame) es critico: un desplazamiento de 1 o 2 nucleotidos cambia completamente la proteina resultante. Por eso las inserciones/deleciones que no son multiplo de 3 (frameshift) suelen ser muy daninas.\n\nEn bioinformatica, el analisis de uso de codones (codon usage) revela preferencias adaptativas de cada organismo y se usa para optimizar la expresion de genes heterologos en biotecnologia."
  },
  {
    id: "codigo_genetico",
    nombre: "Código genético",
    bloque: "bloque5",
    descripcion:
      "Tabla de correspondencia entre los 64 codones posibles y los 20 aminoacidos (mas senales de inicio y parada). Es practicamente universal en todos los seres vivos.\n\nPropiedades fundamentales:\n- Universal: casi identico en bacterias, plantas, animales y hongos (con pocas excepciones en mitocondrias y algunos protistas).\n- Degenerado: multiples codones para un mismo aminoacido (ej: leucina tiene 6 codones).\n- No ambiguo: cada codon codifica un solo aminoacido.\n- Sin solapamiento: cada nucleotido pertenece a un solo codon.\n\nExcepciones notables:\n- Mitocondrias usan UGA para triptofano (no stop).\n- Mycoplasma usa UGA para triptofano.\n- Selenocisteina (aminoacido 21) se incorpora en contextos especificos.\n\nEl codigo genetico fue descifrado en los anos 1960 por Nirenberg, Khorana y Holley (Premio Nobel 1968). Es la tabla Rosetta de la bioinformatica: permite traducir computacionalmente cualquier secuencia de ADN/ARN a proteina."
  },
  {
    id: "degeneracion_codigo",
    nombre: "Degeneración del código genético",
    bloque: "bloque5",
    descripcion:
      "Propiedad del codigo genetico por la cual un mismo aminoacido puede ser codificado por multiples codones diferentes. Esto ocurre porque hay 61 codones para solo 20 aminoacidos.\n\nPatrones de degeneracion:\n- La tercera posicion del codon (wobble position) es la mas variable.\n- Aminoacidos con 6 codones: Leu, Ser, Arg.\n- Aminoacidos con 4 codones: Val, Pro, Thr, Ala, Gly.\n- Aminoacidos con 1 codon: Met (AUG), Trp (UGW).\n\nSignificado biologico:\n- Reduce el impacto de mutaciones: muchas sustituciones en la 3ra posicion son sinonimas (silenciosas).\n- Las mutaciones sinonimas no cambian la proteina pero pueden afectar la eficiencia de traduccion.\n\nEn bioinformatica, la relacion Ka/Ks (mutaciones no sinonimas vs sinonimas) mide la presion selectiva sobre un gen: Ka/Ks > 1 indica seleccion positiva, Ka/Ks < 1 indica seleccion purificadora (gen conservado)."
  },
  {
    id: "codon_inicio",
    nombre: "Codón de inicio (start)",
    bloque: "bloque5",
    descripcion:
      "Codon que senala el punto de inicio de la traduccion en un ARNm. En la gran mayoria de organismos es AUG, que codifica metionina.\n\nEn procariotas:\n- El AUG de inicio es reconocido por la secuencia Shine-Dalgarno (~8 nt upstream) que posiciona al ribosoma.\n- Codones de inicio alternativos: GUG (~14% en E. coli) y UUG (~3%) tambien se usan.\n- La metionina inicial es formil-metionina (fMet).\n\nEn eucariotas:\n- El ribosoma escanea desde el cap 5' hasta encontrar el primer AUG en un contexto favorable (secuencia Kozak).\n- La metionina inicial no esta formilada.\n\nEn bioinformatica, la prediccion de genes requiere identificar correctamente el codon de inicio. En procariotas, se buscan secuencias Shine-Dalgarno upstream de los ATG. Programas como Prodigal y Glimmer estan especializados en esta tarea."
  },
  {
    id: "codon_terminacion",
    nombre: "Codón de terminación (stop)",
    bloque: "bloque5",
    descripcion:
      "Codones que senalan el final de la traduccion. No codifican ningun aminoacido; en su lugar, son reconocidos por factores de liberacion (release factors).\n\nLos tres codones stop del codigo estandar:\n- UAA (ocre): el mas frecuente en E. coli y eucariotas.\n- UAG (ambar): usado en algunos sistemas para incorporar aminoacidos no naturales.\n- UGA (opalo): codon stop mas \"leaky\" (read-through); en algunos contextos codifica selenocisteina.\n\nCuando el ribosoma encuentra un codon stop:\n1. Un factor de liberacion (RF1 o RF2 en procariotas) entra al sitio A.\n2. Se hidroliza el enlace peptidil-ARNt.\n3. El polipeptido se libera.\n4. El ribosoma se disocia.\n\nEn bioinformatica, los marcos de lectura abiertos (ORFs) se definen como secuencias entre un AUG y el primer codon stop en el mismo marco. La busqueda de ORFs es el primer paso en la anotacion de genomas."
  },
  {
    id: "arn_transferencia",
    nombre: "ARN de transferencia (ARNt)",
    bloque: "bloque5",
    descripcion:
      "Pequena molecula de ARN (~76-90 nucleotidos) que funciona como adaptador molecular durante la traduccion, conectando codones del ARNm con sus aminoacidos correspondientes.\n\nEstructura caracteristica en forma de trebol:\n- Brazo aceptor (3'): donde se une el aminoacido (extremo CCA).\n- Brazo del anticodon: contiene el triplete complementario al codon.\n- Brazo D y brazo TψC: contribuyen a la estructura 3D.\n- La estructura real 3D tiene forma de L.\n\nCada aminoacido tiene al menos un ARNt especifico. Las aminoacil-ARNt sintetasas (20 enzimas) cargan cada ARNt con su aminoacido correcto.\n\nHipotesis del wobble (Crick, 1966): la tercera posicion del codon permite un apareamiento flexible con el anticodon, explicando por que no se necesitan 61 ARNt distintos. E. coli tiene 86 genes de ARNt.\n\nEn bioinformatica, tRNAscan-SE es la herramienta estandar para predecir genes de ARNt en genomas."
  },
  {
    id: "anticodon",
    nombre: "Anticodón",
    bloque: "bloque5",
    descripcion:
      "Secuencia de tres nucleotidos en el ARN de transferencia (ARNt) que es complementaria y antiparalela a un codon del ARNm. El apareamiento codon-anticodon asegura que el aminoacido correcto se incorpore durante la traduccion.\n\nMecanismo:\n- El anticodon se lee en direccion 3'→5' del ARNt.\n- Se aparea con el codon del ARNm (5'→3') en el sitio A del ribosoma.\n- El apareamiento debe ser correcto en las dos primeras posiciones.\n- La tercera posicion permite wobble (apareamiento flexible).\n\nEjemplo: el codon AUG (5'-AUG-3') se aparea con el anticodon UAC (3'-UAC-5') del ARNt-Met.\n\nEl wobble en la tercera posicion permite:\n- G del anticodon se aparea con U o C del codon.\n- Inosina (I) del anticodon se aparea con U, C o A del codon.\n\nEste mecanismo reduce el numero total de ARNt necesarios y es una consecuencia directa de la degeneracion del codigo genetico."
  },
  {
    id: "dogma_central",
    nombre: "Dogma central",
    bloque: "bloque5",
    descripcion:
      "Principio formulado por Francis Crick en 1958 (publicado formalmente en 1970) que describe el flujo de informacion genetica en los sistemas biologicos:\n\nADN → ARN → Proteina\n(Replicacion) (Transcripcion) (Traduccion)\n\nTransferencias adicionales conocidas:\n- Transcripcion reversa: ARN → ADN (retrovirus como VIH, retrotransposones).\n- Replicacion del ARN: ARN → ARN (algunos virus de ARN).\n\nTransferencias que NO ocurren:\n- Proteina → ADN\n- Proteina → ARN\n- Proteina → Proteina\n\nEl dogma central es el marco conceptual de toda la biologia molecular y la bioinformatica. Cada paso del dogma tiene herramientas computacionales asociadas: prediccion de genes (ADN), analisis de expresion (ARN), y modelado estructural (Proteina).\n\nCrick aclaro que 'dogma' no significaba certeza absoluta, sino una hipotesis de trabajo que ha resistido decadas de evidencia."
  },

  // === Bloque 6: Tecnologias ADN Recombinante (14) ===
  {
    id: "pcr",
    nombre: "PCR",
    bloque: "bloque6",
    descripcion:
      "Reaccion en Cadena de la Polimerasa (Polymerase Chain Reaction): tecnica inventada por Kary Mullis en 1983 (Premio Nobel 1993) que permite amplificar exponencialmente un fragmento especifico de ADN.\n\nCada ciclo tiene 3 etapas (~2-3 min/ciclo):\n1. Desnaturalizacion (~95°C): separa las dos hebras de ADN.\n2. Alineamiento/Annealing (~55-65°C): los cebadores se unen a sus secuencias complementarias.\n3. Extension (~72°C): la Taq polimerasa (termoestable, de Thermus aquaticus) sintetiza nuevas hebras.\n\nDespues de n ciclos: 2^n copias. 30 ciclos = ~1,000 millones de copias.\n\nVariantes: RT-PCR (detecta ARN), qPCR (cuantitativa, tiempo real), PCR digital (absoluta).\n\nEn bioinformatica, el diseno de cebadores (primers) se realiza con herramientas como Primer3, considerando Tm, contenido GC, dimerización y especificidad contra el genoma diana."
  },
  {
    id: "cebador",
    nombre: "Cebador/Primer",
    bloque: "bloque6",
    descripcion:
      "Oligonucleotido corto (tipicamente 18-25 nucleotidos) de secuencia definida que se aparea por complementariedad con el ADN molde y proporciona un extremo 3'-OH libre para que la ADN polimerasa inicie la sintesis.\n\nParametros de diseno:\n- Longitud: 18-25 nt (compromiso entre especificidad y eficiencia).\n- Contenido GC: 40-60% ideal.\n- Temperatura de melting (Tm): 55-65°C, similar entre forward y reverse.\n- Sin estructuras secundarias (hairpins) ni dimeros primer-primer.\n- Extremo 3' terminando en G o C (\"GC clamp\") para mayor estabilidad.\n\nSe necesitan dos cebadores por PCR: forward (sentido) y reverse (antisentido), flanqueando la region diana.\n\nEn bioinformatica, Primer-BLAST (NCBI) combina el diseno de cebadores con verificacion de especificidad contra genomas completos, evitando amplificaciones inespecificas."
  },
  {
    id: "desnaturalizacion",
    nombre: "Desnaturalización (PCR)",
    bloque: "bloque6",
    descripcion:
      "Primera etapa de cada ciclo de PCR en la que se calienta la mezcla de reaccion a ~94-98°C durante 15-30 segundos para romper los puentes de hidrogeno entre las bases complementarias y separar las dos hebras del ADN.\n\nDetalles:\n- Los pares G-C (3 puentes de H) requieren mas energia que A-T (2 puentes de H).\n- ADN con mayor contenido GC necesita temperaturas mas altas de desnaturalizacion.\n- La desnaturalizacion inicial (antes del primer ciclo) suele ser mas larga (~3-5 min) para asegurar la separacion completa del ADN molde.\n- Temperaturas excesivas o tiempos prolongados danan la Taq polimerasa.\n\nLa Tm (temperatura de melting) de un fragmento de ADN depende de su longitud y contenido GC. Formula basica: Tm ≈ 2(A+T) + 4(G+C) para oligonucleotidos cortos. En bioinformatica, existen calculadoras de Tm mas sofisticadas que consideran la concentracion de sal, mismatches y parametros termodinamicos."
  },
  {
    id: "alineamiento_priming",
    nombre: "Alineamiento/Priming (PCR)",
    bloque: "bloque6",
    descripcion:
      "Segunda etapa del ciclo de PCR en la que se reduce la temperatura (~50-65°C) para permitir que los cebadores se unan (hibriden) con sus secuencias complementarias en el ADN molde.\n\nFactores criticos:\n- La temperatura de annealing debe ser ~5°C por debajo de la Tm del cebador.\n- Temperatura muy baja: los cebadores se unen inespecificamente → bandas multiples.\n- Temperatura muy alta: los cebadores no se unen → no hay amplificacion.\n- El tiempo tipico es 15-30 segundos.\n\nGradiente de temperatura: se puede probar un rango de temperaturas de annealing en paralelo para optimizar la especificidad.\n\nTouchdown PCR: estrategia que comienza con alta temperatura de annealing y la reduce gradualmente, favoreciendo primero la amplificacion especifica."
  },
  {
    id: "extension_pcr",
    nombre: "Extensión (PCR)",
    bloque: "bloque6",
    descripcion:
      "Tercera etapa del ciclo de PCR en la que la ADN polimerasa termoestable sintetiza nuevas hebras complementarias a partir del extremo 3' de cada cebador, usando el ADN molde como guia.\n\nDetalles:\n- Temperatura optima para Taq polimerasa: 72°C.\n- Velocidad de sintesis: ~1,000 nucleotidos/segundo (Taq).\n- Tiempo de extension: ~1 min por cada 1 kb de producto.\n- Para productos largos: se extiende el tiempo de extension.\n\nPolimerasas alternativas:\n- Pfu polimerasa (Pyrococcus furiosus): tiene actividad proofreading (3'→5' exonucleasa), 10x mas fiel que Taq.\n- Phusion: alta fidelidad y velocidad, usada para clonacion.\n\nAl final de la extension, se tiene el doble de copias del fragmento diana. La extension final (despues del ultimo ciclo) suele ser de 5-10 min para completar todos los fragmentos parciales."
  },
  {
    id: "clonacion_molecular",
    nombre: "Clonación molecular",
    bloque: "bloque6",
    descripcion:
      "Conjunto de tecnicas para insertar un fragmento de ADN de interes en un vector replicable y propagarlo en un organismo hospedero (tipicamente E. coli), produciendo millones de copias identicas.\n\nPasos generales:\n1. Obtener el fragmento de ADN (PCR, restriccion, sintesis).\n2. Preparar el vector (cortar con enzimas de restriccion o linearizar).\n3. Ligar el inserto al vector (ADN ligasa).\n4. Transformar celulas competentes con el constructo.\n5. Seleccionar colonias positivas (antibiotico, seleccion azul/blanco).\n6. Verificar el inserto (PCR de colonia, secuenciacion).\n\nAplicaciones: produccion de proteinas recombinantes (ej: insulina), estudios funcionales de genes, construccion de bibliotecas genomicas/de ADNc.\n\nEn la era moderna, la clonacion se complementa con CRISPR-Cas9 para edicion precisa y con sintesis de genes para construir secuencias de novo. En bioinformatica, herramientas como SnapGene y Benchling facilitan el diseno in silico de estrategias de clonacion."
  },
  {
    id: "vector_clonacion",
    nombre: "Vector de clonación",
    bloque: "bloque6",
    descripcion:
      "Molecula de ADN capaz de replicarse autonomamente en una celula hospedero y que permite transportar un fragmento de ADN foraneo (inserto).\n\nTipos principales:\n- Plasmidos: ADN circular, 2-15 kb de inserto. Los mas usados (pUC, pBR322, pET).\n- Bacteriofagos (λ): 15-25 kb de inserto. Utiles para bibliotecas.\n- Cosmidos: hibrido plasmido/fago, 35-45 kb de inserto.\n- BACs (Cromosomas Artificiales Bacterianos): 100-300 kb. Usados en proyectos genoma.\n- YACs (Cromosomas Artificiales de Levadura): 200-2000 kb. Para fragmentos muy grandes.\n\nComponentes esenciales de un vector:\n- Origen de replicacion (ori): permite la replicacion autonoma.\n- Gen de resistencia a antibiotico: seleccion de transformantes.\n- Sitio de clonacion multiple (MCS/polylinker): sitios de restriccion unicos para insertar ADN.\n\nEn bioinformatica, las secuencias de vectores se registran en bases de datos como Addgene para compartir constructos entre laboratorios."
  },
  {
    id: "biblioteca_clones",
    nombre: "Biblioteca de clones",
    bloque: "bloque6",
    descripcion:
      "Coleccion de clones (celulas hospederas con vectores recombinantes) que en conjunto representan fragmentos de un genoma o transcriptoma completo.\n\nTipos:\n- Biblioteca genomica: fragmentos aleatorios del ADN total. Cubre todo el genoma incluyendo intrones, regiones regulatorias e intergenicas.\n- Biblioteca de ADNc: generada a partir de ARNm (via transcripcion reversa). Solo contiene secuencias expresadas (exones).\n- Biblioteca metagenómica: ADN total de una comunidad microbiana (sin cultivar).\n\nCobertura: el numero de clones necesario depende del tamano del genoma y del inserto. Para 99% de probabilidad de tener cualquier secuencia: N = ln(1-0.99)/ln(1-I/G), donde I = tamano inserto, G = tamano genoma.\n\nHistoricamente, las bibliotecas de BAC fueron fundamentales para el Proyecto Genoma Humano (estrategia de shotgun jerarquico). Hoy, la secuenciacion de nueva generacion ha reemplazado muchos usos de las bibliotecas fisicas."
  },
  {
    id: "enzima_restriccion",
    nombre: "Enzima de restricción",
    bloque: "bloque6",
    descripcion:
      "Endonucleasa que reconoce y corta ADN en secuencias especificas (sitios de restriccion). Son parte del sistema inmune bacteriano contra fagos (sistema restriccion-modificacion).\n\nTipos:\n- Tipo II (las mas usadas en laboratorio): reconocen secuencias palindromicas de 4-8 pb y cortan en el sitio o cerca de el. Ej: EcoRI (G|AATTC), BamHI (G|GATCC), HindIII (A|AGCTT).\n- Tipo I y III: cortan lejos del sitio de reconocimiento.\n\nEjemplos con sus cortes:\n- EcoRI: genera extremos pegajosos (5' overhang de 4 nt).\n- SmaI: genera extremos romos (corta en el centro de CCCGGG).\n- Cada enzima tiene temperatura y buffer optimos.\n\nSe han catalogado >4,000 enzimas de restriccion. NEB (New England Biolabs) y herramientas como NEBcutter permiten buscar sitios de restriccion en cualquier secuencia in silico.\n\nSistema CRISPR-Cas: es conceptualmente similar pero programable — puede dirigirse a cualquier secuencia, no solo palindromas."
  },
  {
    id: "sitio_reconocimiento",
    nombre: "Sitio de reconocimiento",
    bloque: "bloque6",
    descripcion:
      "Secuencia especifica de ADN (tipicamente 4-8 pares de bases) donde una enzima de restriccion se une y corta la doble hebra.\n\nCaracteristicas:\n- La mayoria son palindromas (se leen igual en ambas hebras en direccion 5'→3').\n  Ejemplo: EcoRI reconoce 5'-GAATTC-3' / 3'-CTTAAG-5'.\n- Frecuencia esperada: un sitio de 6 bp aparece en promedio cada 4^6 = 4,096 bp.\n- Un sitio de 4 bp aparece cada 4^4 = 256 bp (cortes mas frecuentes).\n- Un sitio de 8 bp aparece cada 4^8 = 65,536 bp (cortes raros, utiles para mapas fisicos).\n\nSitios degenerados: algunas enzimas reconocen secuencias con posiciones ambiguas. Ej: BstEII reconoce G|GTNACC donde N es cualquier base.\n\nEn bioinformatica, la busqueda de sitios de restriccion en secuencias es rutinaria para planificar experimentos de clonacion, mapeo y fingerprinting genomico."
  },
  {
    id: "extremos_romos",
    nombre: "Extremos romos (blunt ends)",
    bloque: "bloque6",
    descripcion:
      "Tipo de extremo generado cuando una enzima de restriccion corta ambas hebras de ADN exactamente en la misma posicion, sin dejar salientes monocatenarios.\n\nEjemplos de enzimas que generan extremos romos:\n- SmaI: CCC|GGG\n- EcoRV: GAT|ATC\n- HaeIII: GG|CC\n\nCaracteristicas:\n- Cualquier extremo romo puede ligarse con cualquier otro extremo romo (no requiere compatibilidad de secuencia).\n- La eficiencia de ligacion es menor que con extremos pegajosos (no hay hibridacion previa que mantenga los fragmentos juntos).\n- Se pueden usar concentraciones mas altas de ADN ligasa o incubaciones mas largas.\n\nAplicaciones: clonacion direccional independiente de secuencia, ligacion de productos de PCR (la Taq polimerasa deja un overhang de A, pero la Pfu deja extremos romos)."
  },
  {
    id: "extremos_pegajosos",
    nombre: "Extremos pegajosos (sticky ends)",
    bloque: "bloque6",
    descripcion:
      "Tipo de extremo generado cuando una enzima de restriccion corta las dos hebras de ADN en posiciones escalonadas, dejando salientes monocatenarios de 1-4 nucleotidos.\n\nTipos:\n- 5' overhang (saliente en 5'): EcoRI genera 5'-G...AATTC-3'. El saliente AATT queda en 5'.\n- 3' overhang (saliente en 3'): KpnI genera 5'-GGTAC...C-3'. El saliente CATG queda en 3'.\n\nVentajas para clonacion:\n- Los salientes complementarios hibridan espontaneamente por puentes de hidrogeno.\n- Esto mantiene los fragmentos unidos temporalmente, facilitando la ligacion.\n- Permiten clonacion direccional: si se usan dos enzimas diferentes, el inserto solo puede entrar en una orientacion.\n\nEnzimas con extremos compatibles: diferentes enzimas pueden generar salientes identicos (ej: BamHI y BglII), permitiendo ligar fragmentos cortados con enzimas distintas (aunque el sitio recombinante ya no es cortable por ninguna de las dos)."
  },
  {
    id: "hibridacion",
    nombre: "Hibridación",
    bloque: "bloque6",
    descripcion:
      "Proceso por el cual dos cadenas complementarias de acidos nucleicos (ADN-ADN, ADN-ARN, o ARN-ARN) se unen formando una molecula de doble hebra mediante puentes de hidrogeno entre bases complementarias.\n\nAplicaciones fundamentales:\n- Southern blot: detecta secuencias especificas de ADN con sondas marcadas.\n- Northern blot: detecta ARN especifico.\n- FISH (Hibridacion in situ fluorescente): localiza secuencias en cromosomas.\n- Microarreglos: miles de hibridaciones simultaneas en un chip.\n- PCR: los cebadores hibridan con el molde.\n\nFactores que afectan la hibridacion:\n- Stringency (rigor): temperatura alta y baja concentracion de sal favorecen apareamiento especifico.\n- Longitud de la sonda: mayor longitud → mayor estabilidad del hibrido.\n- Mismatches: reducen la estabilidad y la Tm del hibrido.\n\nEn bioinformatica, los algoritmos de alineamiento (BLAST, Smith-Waterman) simulan computacionalmente el proceso de hibridacion para encontrar secuencias similares en bases de datos."
  },
  {
    id: "ligacion",
    nombre: "Ligación",
    bloque: "bloque6",
    descripcion:
      "Formacion de un enlace covalente fosfodiester entre el extremo 3'-OH de un fragmento de ADN y el extremo 5'-fosfato del siguiente, sellando el esqueleto azucar-fosfato. Catalizada por la enzima ADN ligasa.\n\nLigasa mas usada: T4 DNA Ligase (del bacteriofago T4).\n- Funciona a 16°C (para extremos pegajosos) o a temperatura ambiente (para extremos romos).\n- Requiere ATP como cofactor.\n- Puede ligar extremos romos y pegajosos.\n\nConsideraciones practicas:\n- Relacion molar inserto:vector optima: 3:1 a 5:1 para extremos pegajosos.\n- Defosforilar el vector (con fosfatasa) previene la re-circularizacion sin inserto.\n- La eficiencia de ligacion de extremos romos es ~100x menor que con extremos pegajosos.\n\nAlternativas modernas: Gibson Assembly, Golden Gate, In-Fusion permiten ligar multiples fragmentos simultaneamente sin necesidad de sitios de restriccion compatibles."
  },

  // === Bloque 7: Tecnicas de Analisis (3) ===
  {
    id: "electroforesis_gel",
    nombre: "Electroforesis en gel",
    bloque: "bloque7",
    descripcion:
      "Tecnica de separacion de moleculas basada en su migracion a traves de una matriz de gel bajo la influencia de un campo electrico. El ADN (cargado negativamente por los fosfatos) migra hacia el polo positivo.\n\nTipos de gel:\n- Agarosa (0.5-2%): separa ADN de 100 bp a 25 kb. Resolucion baja pero rapido y economico.\n- Poliacrilamida (PAGE): separa ADN de 5-500 bp con resolucion de 1 bp. Usado en secuenciacion clasica.\n- SDS-PAGE: separa proteinas por tamano (desnaturalizadas).\n\nVisualizacion:\n- Bromuro de etidio (intercalante, fluorescente bajo UV) — toxicidad, siendo reemplazado.\n- SYBR Safe, GelRed: alternativas menos toxicas.\n- Marcador de peso molecular (ladder): fragmentos de tamano conocido para estimar el tamano de las bandas.\n\nAplicaciones: verificar productos de PCR, evaluar calidad de ADN/ARN extraido, analizar digestiones con enzimas de restriccion. En bioinformatica, la electroforesis in silico simula patrones de bandas a partir de secuencias y sitios de restriccion conocidos."
  },
  {
    id: "sonda",
    nombre: "Sonda (probe)",
    bloque: "bloque7",
    descripcion:
      "Fragmento de acido nucleico (ADN o ARN) de secuencia conocida, marcado con un reporter detectable, que se usa para identificar secuencias complementarias mediante hibridacion.\n\nTipos de marcaje:\n- Radioactivo: P-32, S-35. Alta sensibilidad, pero requiere precauciones de seguridad.\n- Fluorescente: FITC, Cy3, Cy5, Alexa Fluor. Permite multiplexing (multiples colores).\n- Enzimatico: biotina-estreptavidina, digoxigenina. Deteccion colorimetrica.\n\nTipos de sondas:\n- Oligonucleotidos sinteticos: 15-70 nt, diseno preciso.\n- Fragmentos de ADN clonados: 100 bp - varios kb.\n- ARN transcritos in vitro: ribosondas, alta afinidad.\n\nAplicaciones:\n- Southern/Northern blot: deteccion especifica.\n- FISH: localizacion cromosomica.\n- Microarreglos: miles de sondas en un chip.\n\nEn bioinformatica, el diseno de sondas requiere verificar unicidad en el genoma diana (BLAST) y optimizar Tm y especificidad, evitando homologia con secuencias no deseadas."
  },
  {
    id: "microarreglo",
    nombre: "Microarreglo / DNA array",
    bloque: "bloque7",
    descripcion:
      "Tecnologia de alto rendimiento que permite medir simultaneamente la expresion de miles de genes (o detectar variantes genomicas) en un solo experimento.\n\nPrincipio: miles de sondas de ADN de secuencia conocida estan fijadas en posiciones definidas sobre un chip de vidrio o silicio. Se hibrida ARNm (convertido a ADNc marcado con fluorescencia) de la muestra. La intensidad de fluorescencia en cada punto indica el nivel de expresion del gen correspondiente.\n\nTipos:\n- Expresion: mide niveles de ARNm (ej: Affymetrix GeneChip).\n- SNP array: detecta polimorfismos de un solo nucleotido (ej: Illumina Infinium).\n- CGH (Hibridacion Genomica Comparativa): detecta ganancias/perdidas de ADN.\n- ChIP-chip: identifica sitios de union de proteinas al ADN.\n\nAunque RNA-seq ha superado a los microarreglos para muchas aplicaciones, los SNP arrays siguen siendo muy usados en estudios GWAS y en genotipado clinico.\n\nEn bioinformatica, el analisis de microarreglos requiere normalizacion, control de calidad y pruebas estadisticas para identificar genes diferencialmente expresados."
  },

  // === Bloque 8: Variacion y Evolucion (7) ===
  {
    id: "variacion_genetica",
    nombre: "Variación genética intraespecífica",
    bloque: "bloque8",
    descripcion:
      "Diferencias en la secuencia de ADN entre individuos de una misma especie. Es la materia prima sobre la que actua la seleccion natural y otros mecanismos evolutivos.\n\nTipos principales:\n- SNPs (Single Nucleotide Polymorphisms): variaciones de una sola base. Son las mas comunes (~1 SNP cada 1,000 bp en humanos).\n- Indels: inserciones y deleciones pequenas.\n- CNVs (Copy Number Variations): duplicaciones o deleciones de segmentos grandes.\n- Inversiones y translocaciones: reordenamientos estructurales.\n- Microsatelites (STRs): repeticiones cortas en tandem de longitud variable.\n\nFuentes de variacion:\n- Mutacion (fuente primaria).\n- Recombinacion (mezcla variantes existentes).\n- Transferencia horizontal de genes (en bacterias).\n\nEn E. coli, la variacion entre cepas es enorme: el pan-genoma de E. coli contiene >16,000 genes, mientras que cada cepa tiene solo ~4,000-5,500. En bioinformatica, variant calling (GATK, FreeBayes) y GWAS identifican variantes asociadas a fenotipos."
  },
  {
    id: "genoma_referencia",
    nombre: "Genoma de referencia",
    bloque: "bloque8",
    descripcion:
      "Secuencia de ADN completa representativa de una especie, utilizada como patron de comparacion para analisis genomicos. No es el genoma de un individuo \"normal\" sino una construccion consenso.\n\nEjemplos importantes:\n- E. coli K-12 MG1655 (NC_000913.3): ~4.64 Mb, completado en 1997 por Blattner et al.\n- Humano GRCh38/hg38: ~3.1 Gb, version actual del Human Genome Project.\n- Saccharomyces cerevisiae S288C: ~12 Mb, primer eucariota secuenciado (1996).\n\nUsos del genoma de referencia:\n- Mapeo de lecturas de secuenciacion (BWA, Bowtie2).\n- Deteccion de variantes (comparacion muestra vs referencia).\n- Anotacion: identificacion de genes, elementos regulatorios, repeticiones.\n- Base para comparaciones entre cepas y especies.\n\nLimitaciones: un solo genoma no captura toda la diversidad de una especie. Por eso surgen los pan-genomas (union de todos los genes de una especie) y los genomas de referencia pangenomicos (T2T Consortium)."
  },
  {
    id: "conservacion_genetica",
    nombre: "Conservación genética",
    bloque: "bloque8",
    descripcion:
      "Mantenimiento de secuencias de ADN o proteinas con alta similitud entre especies evolutivamente distantes, indicando que la seleccion purificadora ha preservado su funcion durante millones de anos.\n\nNiveles de conservacion:\n- Genes housekeeping: esenciales para la vida, altamente conservados (ej: ARNr 16S, histonas, enzimas del metabolismo central).\n- Dominios proteicos: modulos funcionales conservados que se combinan de diferentes formas.\n- Elementos ultraconservados: regiones con >200 bp identicas entre humano, raton y rata. Muchos son reguladores no codificantes.\n\nMedicion:\n- Identidad de secuencia (% de posiciones identicas en un alineamiento).\n- Ka/Ks ratio: compara mutaciones no-sinonimas vs sinonimas. Ka/Ks << 1 indica fuerte conservacion.\n- Puntuaciones PhyloP y PhastCons: miden conservacion multi-especie.\n\nEn bioinformatica, la conservacion es clave para: predecir funcion de genes desconocidos (ortologos en otras especies ya anotados), identificar regiones regulatorias, y priorizar variantes clinicamente relevantes (las que caen en regiones conservadas son mas probablemente patogenicas)."
  },
  {
    id: "evolucion",
    nombre: "Evolución",
    bloque: "bloque8",
    descripcion:
      "Proceso de cambio en las frecuencias alelicas de una poblacion a traves de generaciones sucesivas, resultando en la diversificacion de los seres vivos.\n\nMecanismos principales:\n- Seleccion natural: variantes mas aptas se reproducen mas.\n- Deriva genetica: cambios aleatorios en poblaciones pequenas.\n- Flujo genico: migracion introduce nuevos alelos.\n- Mutacion: genera nuevas variantes.\n\nEvidencias moleculares:\n- Universalidad del codigo genetico → ancestro comun.\n- Genes homologos entre especies muy distantes.\n- El reloj molecular: la tasa de mutacion neutral permite estimar tiempos de divergencia.\n\nEn bacterias, la evolucion es particularmente rapida por:\n- Generaciones cortas (~20 min en E. coli).\n- Transferencia horizontal de genes (conjugacion, transformacion, transduccion).\n- Grandes tamanos poblacionales.\n\nEn bioinformatica, la filogenómica reconstruye la historia evolutiva comparando genomas completos. Herramientas como RAxML, IQ-TREE y BEAST construyen arboles filogeneticos con modelos estadisticos sofisticados."
  },
  {
    id: "seleccion_natural",
    nombre: "Selección natural",
    bloque: "bloque8",
    descripcion:
      "Mecanismo evolutivo propuesto por Charles Darwin (1859) por el cual los individuos con variantes geneticas que confieren mayor aptitud reproductiva (fitness) dejan mas descendencia, aumentando la frecuencia de esas variantes en la poblacion.\n\nTipos:\n- Seleccion positiva (direccional): favorece nuevas variantes ventajosas. Ka/Ks > 1.\n- Seleccion purificadora (negativa): elimina variantes daninas. Ka/Ks < 1. La mas comun.\n- Seleccion balanceadora: mantiene multiples alelos (ej: heterocigoto superior).\n- Seleccion diversificadora: favorece fenotipos extremos.\n\nEjemplos clasicos en bacterias:\n- Resistencia a antibioticos: mutaciones en genes diana o adquisicion de genes de resistencia.\n- Seleccion de cepas patogenas con factores de virulencia.\n\nEn bioinformatica, la deteccion de seleccion se realiza mediante:\n- Test Ka/Ks (dN/dS) entre ortologos.\n- Test de McDonald-Kreitman.\n- Analisis de barrido selectivo (selective sweep) en datos poblacionales.\n- iHS y Fst para detectar seleccion reciente en genomas humanos."
  },
  {
    id: "adaptacion",
    nombre: "Adaptación",
    bloque: "bloque8",
    descripcion:
      "Proceso por el cual una poblacion se ajusta geneticamente a su ambiente a traves de la seleccion natural, incrementando la frecuencia de rasgos heredables que mejoran la supervivencia y reproduccion.\n\nNiveles de adaptacion:\n- Molecular: cambios en la secuencia de proteinas que mejoran su funcion en el ambiente (ej: enzimas termófilas en organismos de fuentes termales).\n- Fisiologica: cambios en la regulacion genica (ej: induccion de genes de resistencia a antibioticos).\n- Morfologica: cambios estructurales (ej: diferencias en picos de pinzones de Darwin).\n\nAdaptacion en E. coli:\n- El experimento LTEE (Long-Term Evolution Experiment) de Richard Lenski: 12 poblaciones de E. coli evolucionando desde 1988 (>75,000 generaciones). Se han observado adaptaciones como la capacidad de metabolizar citrato, cambios en tamano celular y tasas de mutacion.\n\nEn bioinformatica, la genomica comparativa entre cepas adaptadas a diferentes nichos revela los genes y mutaciones responsables de la adaptacion, mediante analisis de convergencia evolutiva y seleccion positiva."
  },
  {
    id: "especiacion",
    nombre: "Especiación",
    bloque: "bloque8",
    descripcion:
      "Proceso por el cual una poblacion ancestral diverge en dos o mas especies reproductivamente aisladas. Es el mecanismo que genera la biodiversidad.\n\nModos principales:\n- Alopatrica: separacion geografica impide el flujo genico. La mas comun.\n- Simpatrica: divergencia sin barrera fisica, por seleccion disruptiva o especializacion ecologica.\n- Parapatrica: divergencia en zonas adyacentes con flujo genico limitado.\n\nEn bacterias, el concepto de especie es diferente:\n- No hay reproduccion sexual clasica.\n- Se usan criterios genomicos: ANI (Average Nucleotide Identity) >95% = misma especie.\n- DDH (DNA-DNA Hybridization) >70% fue el estandar clasico.\n- El analisis del gen 16S rRNA >97% de identidad como criterio minimo.\n\nEn bioinformatica, la delimitacion de especies bacterianas se basa en analisis de genoma completo: ANI (pyani, fastANI), AAI (Average Amino acid Identity), y arboles filogeneticos de genes concatenados (core genome MLSA)."
  },

  // === Bloque 9: Bioinformatica (5) ===
  {
    id: "alfabeto_molecular",
    nombre: "Alfabeto molecular",
    bloque: "bloque9",
    descripcion:
      "Representacion de macromoleculas biologicas como cadenas de caracteres de un alfabeto finito, lo que permite su manipulacion computacional.\n\nAlfabetos:\n- ADN: {A, T, G, C} — 4 letras. Con ambiguedades IUPAC: R(AG), Y(CT), M(AC), K(GT), S(GC), W(AT), N(cualquiera).\n- ARN: {A, U, G, C} — 4 letras.\n- Proteinas: {A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y} — 20 letras (codigo de 1 letra IUPAC).\n\nEsta representacion es la piedra angular de la bioinformatica:\n- Permite almacenar genomas como archivos de texto (formato FASTA).\n- Habilita algoritmos de busqueda de patrones (expresiones regulares, motifs).\n- Permite alineamientos de secuencias (comparacion caracter a caracter).\n- Posibilita calculos estadisticos (frecuencia de letras, entropia, composicion).\n\nEl genoma de E. coli K-12 como texto: ~4.6 millones de caracteres, ocupando ~4.6 MB en formato plano. Todo el conocimiento genetico de un organismo codificado en 4 letras."
  },
  {
    id: "genomica_comparativa",
    nombre: "Genómica comparativa",
    bloque: "bloque9",
    descripcion:
      "Subdisciplina de la bioinformatica que compara genomas completos de diferentes organismos para inferir relaciones evolutivas, identificar genes conservados y descubrir diferencias funcionales.\n\nAnalisis tipicos:\n- Identificacion de ortologos y paralogos (genes homologos entre/dentro de especies).\n- Sintenia: conservacion del orden de genes entre genomas.\n- Pan-genoma: core genome (genes compartidos por todas las cepas) + accessory genome (genes de algunas cepas) + unique genes.\n- Islas genomicas: regiones presentes en una cepa pero ausentes en otras (frecuentemente adquiridas por transferencia horizontal).\n- ANI (Average Nucleotide Identity): medida de similitud global entre dos genomas.\n\nHerramientas clave:\n- OrthoFinder, Roary, Panaroo: analisis de ortologos y pan-genoma.\n- Mauve, MUMmer, Minimap2: alineamiento de genomas completos.\n- Circos, Artemis: visualizacion de comparaciones.\n\nEn este proyecto (GenomeHub), la comparacion genomica entre cepas de E. coli permite entender que genes comparten, cuales son unicos, y como han divergido evolutivamente."
  },
  {
    id: "alineamiento_secuencias",
    nombre: "Alineamiento de secuencias",
    bloque: "bloque9",
    descripcion:
      "Metodo computacional fundamental de la bioinformatica que consiste en disponer dos o mas secuencias (de ADN, ARN o proteina) una sobre otra, insertando gaps cuando es necesario, para maximizar las coincidencias y revelar regiones de similitud.\n\nTipos:\n- Pairwise (2 secuencias):\n  · Global (Needleman-Wunsch, 1970): alinea las secuencias completas extremo a extremo.\n  · Local (Smith-Waterman, 1981): encuentra la region de mayor similitud local.\n- Multiple (MSA, >2 secuencias):\n  · ClustalW/Omega, MUSCLE, MAFFT, T-Coffee.\n  · Esencial para filogenetica y deteccion de motifs conservados.\n\nComponentes del scoring:\n- Match: puntuacion positiva por coincidencia.\n- Mismatch: penalizacion por sustitucion.\n- Gap penalty: penalizacion por insercion/delecion (apertura + extension).\n- Matrices de sustitucion: BLOSUM62 (proteinas), PAM (proteinas), NUC.4.4 (ADN).\n\nEl alineamiento es la operacion mas basica y frecuente en bioinformatica: se usa para buscar genes homologos, predecir funcion, detectar mutaciones, construir arboles filogeneticos y validar ensamblajes genomicos."
  },
  {
    id: "blast",
    nombre: "BLAST",
    bloque: "bloque9",
    descripcion:
      "Basic Local Alignment Search Tool: familia de algoritmos desarrollados por Altschul et al. (1990) en el NCBI para buscar rapidamente secuencias similares en bases de datos masivas.\n\nVariantes:\n- BLASTn: nucleotido vs nucleotido.\n- BLASTp: proteina vs proteina.\n- BLASTx: nucleotido (traducido en 6 marcos) vs proteina.\n- tBLASTn: proteina vs nucleotido (traducido en 6 marcos).\n- tBLASTx: nucleotido traducido vs nucleotido traducido.\n\nAlgoritmo (simplificado):\n1. Genera palabras cortas (words) de la query.\n2. Busca coincidencias exactas (seeds) en la base de datos.\n3. Extiende los seeds en ambas direcciones (HSPs).\n4. Evalua significancia estadistica (E-value).\n\nE-value: numero esperado de alineamientos con score igual o mejor por azar. E-value < 1e-5 se considera significativo.\n\nBLAST procesa miles de queries por segundo contra bases de datos de millones de secuencias. Es la herramienta mas citada de la bioinformatica y la puerta de entrada al analisis de cualquier secuencia desconocida."
  },
  {
    id: "bioinformatica",
    nombre: "Bioinformática",
    bloque: "bloque9",
    descripcion:
      "Disciplina interdisciplinaria que aplica matematicas, estadistica y ciencias de la computacion al analisis de datos biologicos, especialmente secuencias de ADN, ARN y proteinas.\n\nAreas principales:\n- Genomica: ensamblaje, anotacion y comparacion de genomas.\n- Transcriptomica: analisis de expresion genica (RNA-seq, microarreglos).\n- Proteomica: identificacion y cuantificacion de proteinas (espectrometria de masas).\n- Filogenetica: reconstruccion de relaciones evolutivas.\n- Biologia estructural: prediccion de estructura 3D de proteinas (AlphaFold).\n- Metagenómica: analisis de comunidades microbianas sin cultivo.\n- Farmacogenómica: variacion genetica y respuesta a farmacos.\n\nHitos historicos:\n- 1965: Atlas of Protein Sequence (Margaret Dayhoff) → primera base de datos.\n- 1970: Needleman-Wunsch → primer algoritmo de alineamiento.\n- 1990: BLAST → busqueda rapida en bases de datos.\n- 2001: Genoma Humano completado.\n- 2022: AlphaFold2 predice estructura de ~200 millones de proteinas.\n\nBases de datos fundamentales: GenBank (NCBI), UniProt (proteinas), PDB (estructuras 3D), KEGG (vias metabolicas), Ensembl (genomas anotados)."
  }
];

// Links: conexiones lógicas (~200) con verbo en `relacion`
const LINKS = [
  // --- Bloque 1 interno y con macromoléculas ---
  { source: "celula", target: "teoria_celular", relacion: "se fundamenta en" },
  { source: "celula", target: "ciclo_vida_celular", relacion: "atraviesa" },
  { source: "celula", target: "via_metabolica", relacion: "ejecuta" },
  { source: "ciclo_vida_celular", target: "via_metabolica", relacion: "es regulado por" },
  { source: "ciclo_vida_celular", target: "replicacion_adn", relacion: "requiere" },
  { source: "celula", target: "adn", relacion: "contiene" },
  { source: "celula", target: "arn", relacion: "produce" },
  { source: "celula", target: "proteina", relacion: "utiliza" },
  { source: "via_metabolica", target: "proteina", relacion: "es catalizada por" },

  // --- Bloque 2 interno ---
  { source: "adn", target: "nucleotido", relacion: "se compone de" },
  { source: "arn", target: "nucleotido", relacion: "se compone de" },
  { source: "nucleotido", target: "bases_nitrogenadas", relacion: "incluye" },
  { source: "adn", target: "regla_chargaff", relacion: "obedece a" },
  { source: "regla_chargaff", target: "complementariedad_bases", relacion: "describe" },
  { source: "adn", target: "doble_helice", relacion: "adopta forma de" },
  { source: "doble_helice", target: "complementariedad_bases", relacion: "se estabiliza por" },
  { source: "complementariedad_bases", target: "replicacion_adn", relacion: "permite" },
  { source: "replicacion_adn", target: "mutacion", relacion: "puede introducir" },
  { source: "arn", target: "bases_nitrogenadas", relacion: "utiliza" },
  { source: "proteina", target: "aminoacido", relacion: "se construye a partir de" },

  // --- Bloque 3 interno ---
  { source: "cromosoma", target: "adn", relacion: "organiza" },
  { source: "cromosoma", target: "gen", relacion: "contiene" },
  { source: "gen", target: "adn", relacion: "es segmento de" },
  { source: "herencia_mendeliana", target: "gen", relacion: "describe transmisión de" },
  { source: "herencia_mendeliana", target: "cromosoma", relacion: "se observa en" },
  { source: "ligamiento_genetico", target: "cromosoma", relacion: "depende de su posición en" },
  { source: "mapa_genetico", target: "ligamiento_genetico", relacion: "se infiere a partir de" },
  { source: "mapa_genetico", target: "cromosoma", relacion: "ordena genes en" },
  { source: "un_gen_una_proteina", target: "gen", relacion: "asocia" },
  { source: "un_gen_una_proteina", target: "proteina", relacion: "relaciona con" },
  { source: "mutacion", target: "gen", relacion: "modifica" },
  { source: "mutacion", target: "variacion_genetica", relacion: "genera" },

  // --- Bloque 4 interno ---
  { source: "celula_eucariota", target: "celula", relacion: "es tipo de" },
  { source: "celula_procariota", target: "celula", relacion: "es tipo de" },
  { source: "celula_eucariota", target: "exon", relacion: "presenta genes con" },
  { source: "celula_eucariota", target: "intron", relacion: "presenta genes con" },
  { source: "exon", target: "gen", relacion: "forma parte de" },
  { source: "intron", target: "gen", relacion: "interrumpe a" },
  { source: "celula_procariota", target: "intron", relacion: "carece de" },

  // --- Flujo de información: conexiones Bloque 2–5–4–3 ---
  { source: "replicacion_adn", target: "ciclo_vida_celular", relacion: "sostiene" },
  { source: "replicacion_adn", target: "celula_procariota", relacion: "ocurre en" },
  { source: "replicacion_adn", target: "celula_eucariota", relacion: "ocurre en" },

  // --- Bloque 5 interno (transcripción y procesamiento) ---
  { source: "transcripcion", target: "adn", relacion: "lee" },
  { source: "transcripcion", target: "arn_mensajero", relacion: "genera" },
  { source: "arn_polimerasa", target: "transcripcion", relacion: "cataliza" },
  { source: "arn_polimerasa", target: "nucleotido", relacion: "incorpora" },
  { source: "arn_mensajero", target: "exon", relacion: "conserva" },
  { source: "arn_mensajero", target: "intron", relacion: "excluye tras" },
  { source: "splicing", target: "intron", relacion: "elimina" },
  { source: "splicing", target: "exon", relacion: "une" },
  { source: "splicing", target: "arn_mensajero", relacion: "madura" },
  { source: "dogma_central", target: "transcripcion", relacion: "incluye" },

  // --- Bloque 5 interno (traducción) ---
  { source: "traduccion", target: "arn_mensajero", relacion: "lee" },
  { source: "traduccion", target: "ribosoma", relacion: "es realizada por" },
  { source: "traduccion", target: "proteina", relacion: "sintetiza" },
  { source: "traduccion", target: "polipeptido", relacion: "ensambla" },
  { source: "ribosoma", target: "arn_mensajero", relacion: "se acopla a" },
  { source: "ribosoma", target: "arn_transferencia", relacion: "coordina con" },
  { source: "arn_transferencia", target: "aminoacido", relacion: "transporta" },
  { source: "arn_transferencia", target: "anticodon", relacion: "contiene" },
  { source: "anticodon", target: "codon", relacion: "reconoce" },
  { source: "codon", target: "codigo_genetico", relacion: "se interpreta mediante" },
  { source: "codigo_genetico", target: "aminoacido", relacion: "especifica" },
  { source: "degeneracion_codigo", target: "codigo_genetico", relacion: "caracteriza" },
  { source: "codon_inicio", target: "traduccion", relacion: "inicia" },
  { source: "codon_terminacion", target: "traduccion", relacion: "finaliza" },
  { source: "aminoacido", target: "polipeptido", relacion: "forma" },
  { source: "polipeptido", target: "proteina", relacion: "se pliega en" },
  { source: "dogma_central", target: "traduccion", relacion: "incluye" },
  { source: "dogma_central", target: "proteina", relacion: "culmina en" },

  // --- Bloque 6 interno: PCR y enzimas ---
  { source: "pcr", target: "adn", relacion: "amplifica" },
  { source: "pcr", target: "cebador", relacion: "requiere" },
  { source: "pcr", target: "desnaturalizacion", relacion: "incluye etapa de" },
  { source: "pcr", target: "alineamiento_priming", relacion: "incluye etapa de" },
  { source: "pcr", target: "extension_pcr", relacion: "incluye etapa de" },
  { source: "desnaturalizacion", target: "doble_helice", relacion: "separa" },
  { source: "alineamiento_priming", target: "complementariedad_bases", relacion: "aprovecha" },
  { source: "extension_pcr", target: "replicacion_adn", relacion: "imita" },
  { source: "extension_pcr", target: "nucleotido", relacion: "incorpora" },

  { source: "clonacion_molecular", target: "vector_clonacion", relacion: "usa" },
  { source: "clonacion_molecular", target: "enzima_restriccion", relacion: "aplica" },
  { source: "clonacion_molecular", target: "ligacion", relacion: "requiere" },
  { source: "enzima_restriccion", target: "sitio_reconocimiento", relacion: "reconoce" },
  { source: "enzima_restriccion", target: "extremos_romos", relacion: "puede generar" },
  { source: "enzima_restriccion", target: "extremos_pegajosos", relacion: "puede generar" },
  { source: "extremos_pegajosos", target: "hibridacion", relacion: "facilitan" },
  { source: "hibridacion", target: "ligacion", relacion: "precede a" },
  { source: "ligacion", target: "vector_clonacion", relacion: "sella inserto en" },
  { source: "vector_clonacion", target: "biblioteca_clones", relacion: "genera" },

  // --- Bloque 7 interno: análisis ---
  { source: "electroforesis_gel", target: "pcr", relacion: "se usa para visualizar productos de" },
  { source: "electroforesis_gel", target: "clonacion_molecular", relacion: "verifica" },
  { source: "sonda", target: "hibridacion", relacion: "utiliza" },
  { source: "sonda", target: "microarreglo", relacion: "constituye puntos en" },
  { source: "microarreglo", target: "arn_mensajero", relacion: "mide abundancia de" },
  { source: "microarreglo", target: "dogma_central", relacion: "mapea expresión a lo largo de" },

  // --- Bloque 8 interno: variación y evolución ---
  { source: "variacion_genetica", target: "evolucion", relacion: "alimenta" },
  { source: "variacion_genetica", target: "genoma_referencia", relacion: "se define respecto a" },
  { source: "variacion_genetica", target: "seleccion_natural", relacion: "es filtrada por" },
  { source: "genoma_referencia", target: "conservacion_genetica", relacion: "permite detectar" },
  { source: "conservacion_genetica", target: "evolucion", relacion: "evidencia historia de" },
  { source: "seleccion_natural", target: "adaptacion", relacion: "favorece" },
  { source: "adaptacion", target: "evolucion", relacion: "se manifiesta en" },
  { source: "evolucion", target: "especiacion", relacion: "puede conducir a" },

  // --- Bloque 9 interno: bioinformática ---
  { source: "alfabeto_molecular", target: "adn", relacion: "representa como cadena de" },
  { source: "alfabeto_molecular", target: "proteina", relacion: "representa secuencia de" },
  { source: "alfabeto_molecular", target: "alineamiento_secuencias", relacion: "permite codificar para" },
  { source: "alineamiento_secuencias", target: "variacion_genetica", relacion: "detecta" },
  { source: "alineamiento_secuencias", target: "conservacion_genetica", relacion: "revela" },
  { source: "blast", target: "alineamiento_secuencias", relacion: "acelera cálculo de" },
  { source: "blast", target: "genomica_comparativa", relacion: "apoya" },
  { source: "genomica_comparativa", target: "evolucion", relacion: "reconstruye historia de" },
  { source: "genomica_comparativa", target: "genoma_referencia", relacion: "usa como marco" },
  { source: "bioinformatica", target: "blast", relacion: "desarrolla y aplica" },
  { source: "bioinformatica", target: "alineamiento_secuencias", relacion: "optimiza métodos de" },
  { source: "bioinformatica", target: "microarreglo", relacion: "analiza datos de" },
  { source: "bioinformatica", target: "pcr", relacion: "diseña cebadores para" },
  { source: "bioinformatica", target: "genomica_comparativa", relacion: "implementa algoritmos de" },
  { source: "bioinformatica", target: "dogma_central", relacion: "modela digitalmente" },
  { source: "bioinformatica", target: "via_metabolica", relacion: "reconstruye redes de" },

  // --- Conexiones cruzadas extra para enriquecer el grafo ---

  // Fundamentos ↔ Organización genómica
  { source: "teoria_celular", target: "herencia_mendeliana", relacion: "proporciona marco para" },
  { source: "celula", target: "cromosoma", relacion: "aloja" },
  { source: "ciclo_vida_celular", target: "mutacion", relacion: "puede acumular" },

  // Macromoléculas ↔ Flujo
  { source: "adn", target: "dogma_central", relacion: "inicia flujo descrito por" },
  { source: "arn", target: "dogma_central", relacion: "intermedia en" },
  { source: "proteina", target: "via_metabolica", relacion: "ejecuta reacciones en" },
  { source: "nucleotido", target: "pcr", relacion: "es insumo para" },
  { source: "bases_nitrogenadas", target: "mutacion", relacion: "cambian en" },
  { source: "complementariedad_bases", target: "hibridacion", relacion: "hace posible" },

  // Células ↔ Flujo
  { source: "celula_eucariota", target: "transcripcion", relacion: "realiza en núcleo" },
  { source: "celula_eucariota", target: "traduccion", relacion: "realiza en cytoplasma" },
  { source: "celula_procariota", target: "transcripcion", relacion: "acopla con traducción" },
  { source: "celula_procariota", target: "traduccion", relacion: "ocurre simultáneamente con" },

  // Flujo ↔ Tecnologías
  { source: "pcr", target: "mutacion", relacion: "puede detectar variantes como" },
  { source: "pcr", target: "variacion_genetica", relacion: "amplifica regiones con" },
  { source: "clonacion_molecular", target: "gen", relacion: "manipula copias de" },
  { source: "biblioteca_clones", target: "genomica_comparativa", relacion: "alimenta estudios de" },
  { source: "electroforesis_gel", target: "mutacion", relacion: "permite observar efectos de" },
  { source: "microarreglo", target: "expresion_genica", relacion: "perfilan niveles de" },

  // Variación/Evolución ↔ Bioinformática
  { source: "variacion_genetica", target: "alineamiento_secuencias", relacion: "se revela mediante" },
  { source: "genoma_referencia", target: "blast", relacion: "sirve como base en" },
  { source: "conservacion_genetica", target: "genomica_comparativa", relacion: "se cuantifica con" },
  { source: "evolucion", target: "bioinformatica", relacion: "se estudia con herramientas de" },
  { source: "seleccion_natural", target: "genomica_comparativa", relacion: "deja huellas en" },
  { source: "adaptacion", target: "genomica_comparativa", relacion: "se infiere mediante" },
  { source: "especiacion", target: "genomica_comparativa", relacion: "se analiza comparando" },

  // Más conexiones densas para enriquecer el grafo

  // Entre flujo y variación
  { source: "mutacion", target: "codigo_genetico", relacion: "altera lectura de" },
  { source: "mutacion", target: "proteina", relacion: "puede cambiar estructura de" },
  { source: "mutacion", target: "splicing", relacion: "puede afectar" },
  { source: "variacion_genetica", target: "proteina", relacion: "produce isoformas de" },

  // Flujo interno extra
  { source: "arn_mensajero", target: "codigo_genetico", relacion: "codifica usando" },
  { source: "arn_mensajero", target: "codon", relacion: "se organiza en" },
  { source: "codon", target: "aminoacido", relacion: "especifica" },
  { source: "arn_transferencia", target: "codigo_genetico", relacion: "traduce reglas de" },
  { source: "dogma_central", target: "replicacion_adn", relacion: "se complementa con" },

  // Tecnologías ↔ Bioinformática
  { source: "pcr", target: "bioinformatica", relacion: "se optimiza mediante" },
  { source: "enzima_restriccion", target: "bioinformatica", relacion: "se diseña usando" },
  { source: "microarreglo", target: "bioinformatica", relacion: "genera matrices analizadas por" },
  { source: "electroforesis_gel", target: "bioinformatica", relacion: "puede cuantificarse con" },

  // Fundamentos ↔ Bioinformática
  { source: "teoria_celular", target: "bioinformatica", relacion: "inspira modelos multiescala en" },
  { source: "via_metabolica", target: "genomica_comparativa", relacion: "se compara entre especies mediante" },
  { source: "via_metabolica", target: "bioinformatica", relacion: "se representa como red en" },

  // Organización genómica ↔ Bioinformática
  { source: "cromosoma", target: "genomica_comparativa", relacion: "es unidad de comparación en" },
  { source: "mapa_genetico", target: "bioinformatica", relacion: "se construye con algoritmos de" },
  { source: "ligamiento_genetico", target: "variacion_genetica", relacion: "depende de patrones de" },

  // Bioinformática interna extra
  { source: "alfabeto_molecular", target: "blast", relacion: "alimenta bases de datos para" },
  { source: "alfabeto_molecular", target: "genomica_comparativa", relacion: "permite codificar genomas en" },
  { source: "bioinformatica", target: "variacion_genetica", relacion: "modela el impacto de" },
  { source: "bioinformatica", target: "seleccion_natural", relacion: "detecta firmas de" }
];

// Exportar objeto global para app.js
window.BIOINFO_GRAPH = {
  bloques: BLOQUES,
  nodes: NODOS,
  links: LINKS
};


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

// Nodos: 69 conceptos con bloque y descripción
const NODOS = [
  // === Bloque 1: Fundamentos Celulares (4) ===
  {
    id: "celula",
    nombre: "Célula",
    bloque: "bloque1",
    descripcion:
      "Unidad estructural y funcional básica de los organismos, compuesta por componentes que permiten metabolismo, replicación y respuesta al entorno."
  },
  {
    id: "teoria_celular",
    nombre: "Teoría celular",
    bloque: "bloque1",
    descripcion:
      "Postulado biológico que establece que los organismos están formados por células y que la célula es la unidad fundamental de vida."
  },
  {
    id: "ciclo_vida_celular",
    nombre: "Ciclo de vida celular",
    bloque: "bloque1",
    descripcion:
      "Secuencia de estados por los que pasa una célula (crecimiento, replicación/división y muerte), regulada por procesos bioquímicos."
  },
  {
    id: "via_metabolica",
    nombre: "Vía metabólica / Pathway",
    bloque: "bloque1",
    descripcion:
      "Red organizada de reacciones químicas catalizadas que sintetizan, degradan o señalizan procesos dentro de la célula."
  },

  // === Bloque 2: Macromoléculas Informacionales (9) ===
  {
    id: "adn",
    nombre: "ADN",
    bloque: "bloque2",
    descripcion:
      "Polímero de nucleótidos con bases A, T, G y C; normalmente bicatenario y químicamente estable; almacena información hereditaria."
  },
  {
    id: "arn",
    nombre: "ARN",
    bloque: "bloque2",
    descripcion:
      "Polímero de nucleótidos con bases A, U, G y C; usualmente monocatenario y más reactivo; participa en transferencia y procesamiento de información genética y en funciones catalíticas/regulatorias."
  },
  {
    id: "proteina",
    nombre: "Proteína",
    bloque: "bloque2",
    descripcion:
      "Macromolécula formada por una o más cadenas polipeptídicas de aminoácidos; ejecuta funciones celulares (catálisis, estructura, transporte, señalización, regulación)."
  },
  {
    id: "nucleotido",
    nombre: "Nucleótido",
    bloque: "bloque2",
    descripcion:
      "Unidad básica de ADN/ARN compuesta por una base nitrogenada, un azúcar y un grupo fosfato."
  },
  {
    id: "bases_nitrogenadas",
    nombre: "Bases nitrogenadas",
    bloque: "bloque2",
    descripcion:
      "A (adenina), T (timina, solo ADN), U (uracilo, solo ARN), G (guanina), C (citosina): bases que codifican la información en ácidos nucleicos."
  },
  {
    id: "regla_chargaff",
    nombre: "Regla de Chargaff",
    bloque: "bloque2",
    descripcion:
      "Observación de que en ADN bicatenario la cantidad de A≈T y G≈C, consistente con apareamiento complementario."
  },
  {
    id: "doble_helice",
    nombre: "Doble hélice",
    bloque: "bloque2",
    descripcion:
      "Estructura del ADN formada por dos hebras antiparalelas enrolladas, unidas por puentes de hidrógeno entre pares de bases."
  },
  {
    id: "complementariedad_bases",
    nombre: "Complementariedad de bases",
    bloque: "bloque2",
    descripcion:
      "Emparejamiento específico A–T (o A–U en ARN) y C–G; permite copia fiel en replicación/transcripción e hibridación."
  },
  {
    id: "replicacion_adn",
    nombre: "Replicación del ADN",
    bloque: "bloque2",
    descripcion:
      "Proceso de copia del ADN donde cada hebra sirve como molde para sintetizar una hebra complementaria."
  },

  // === Bloque 3: Organización Genómica (7) ===
  {
    id: "cromosoma",
    nombre: "Cromosoma",
    bloque: "bloque3",
    descripcion:
      "Estructura que organiza y empaqueta ADN (y proteínas asociadas) y contiene genes; su número varía entre especies."
  },
  {
    id: "gen",
    nombre: "Gen",
    bloque: "bloque3",
    descripcion:
      "Segmento de ADN que contiene la información necesaria para producir un producto funcional (ARN o proteína) bajo un contexto regulatorio."
  },
  {
    id: "herencia_mendeliana",
    nombre: "Herencia mendeliana",
    bloque: "bloque3",
    descripcion:
      "Transmisión de rasgos mediante unidades discretas (genes) que segregan y se combinan según reglas observables en la descendencia."
  },
  {
    id: "mutacion",
    nombre: "Mutación",
    bloque: "bloque3",
    descripcion:
      "Cambio en la secuencia de ADN (p. ej., sustitución de base) que puede alterar un rasgo o función biológica."
  },
  {
    id: "ligamiento_genetico",
    nombre: "Ligamiento genético",
    bloque: "bloque3",
    descripcion:
      "Tendencia de genes cercanos en un cromosoma a heredarse juntos, debido a menor probabilidad de recombinación entre ellos."
  },
  {
    id: "mapa_genetico",
    nombre: "Mapa genético",
    bloque: "bloque3",
    descripcion:
      "Ordenamiento relativo de genes en un cromosoma estimado a partir de frecuencias de recombinación (distancia genética)."
  },
  {
    id: "un_gen_una_proteina",
    nombre: "“Un gen–una proteína”",
    bloque: "bloque3",
    descripcion:
      "Hipótesis histórica que asocia cada gen con la producción de una proteína; hoy se reconoce que un gen puede originar múltiples productos."
  },

  // === Bloque 4: Células y Organización (4) ===
  {
    id: "celula_eucariota",
    nombre: "Célula eucariota",
    bloque: "bloque4",
    descripcion:
      "Célula con ADN encapsulado en un núcleo; en general, genes con exones e intrones y procesamiento de ARN."
  },
  {
    id: "celula_procariota",
    nombre: "Célula procariota",
    bloque: "bloque4",
    descripcion:
      "Célula sin núcleo; genes típicamente continuos; transcripción y traducción ocurren en el mismo compartimento."
  },
  {
    id: "exon",
    nombre: "Exón",
    bloque: "bloque4",
    descripcion:
      "Segmento de un gen eucariota que permanece en el ARNm maduro y contribuye a la secuencia codificante."
  },
  {
    id: "intron",
    nombre: "Intrón",
    bloque: "bloque4",
    descripcion:
      "Segmento interveniente transcrito en ARN pero eliminado durante el splicing antes de formar el ARNm maduro."
  },

  // === Bloque 5: Flujo de Información Genética (16) ===
  {
    id: "transcripcion",
    nombre: "Transcripción",
    bloque: "bloque5",
    descripcion:
      "Síntesis de ARN a partir de un molde de ADN mediante una enzima; copia información de un gen a una molécula de ARN."
  },
  {
    id: "arn_polimerasa",
    nombre: "ARN polimerasa",
    bloque: "bloque5",
    descripcion:
      "Complejo enzimático que cataliza la transcripción, agregando ribonucleótidos complementarios al molde de ADN."
  },
  {
    id: "arn_mensajero",
    nombre: "ARN mensajero (ARNm)",
    bloque: "bloque5",
    descripcion:
      "ARN que porta la información codificante de un gen hacia el ribosoma para síntesis proteica."
  },
  {
    id: "splicing",
    nombre: "Splicing (empalme)",
    bloque: "bloque5",
    descripcion:
      "Proceso eucariota de eliminación de intrones y unión de exones para generar ARNm funcional."
  },
  {
    id: "traduccion",
    nombre: "Traducción",
    bloque: "bloque5",
    descripcion:
      "Proceso por el cual el ribosoma lee codones del ARNm y ensambla una proteína incorporando aminoácidos en orden."
  },
  {
    id: "ribosoma",
    nombre: "Ribosoma",
    bloque: "bloque5",
    descripcion:
      "Complejo ribonucleoproteico (ARN + proteínas) que cataliza la traducción y ensambla polipéptidos."
  },
  {
    id: "aminoacido",
    nombre: "Aminoácido",
    bloque: "bloque5",
    descripcion:
      "Monómero de las proteínas; existen 20 tipos estándar codificados por el código genético."
  },
  {
    id: "polipeptido",
    nombre: "Polipéptido",
    bloque: "bloque5",
    descripcion:
      "Cadena lineal de aminoácidos unida por enlaces peptídicos; puede plegarse para formar una proteína funcional."
  },
  {
    id: "codon",
    nombre: "Codón",
    bloque: "bloque5",
    descripcion:
      "Triplete de nucleótidos en ARNm que especifica un aminoácido o una señal de inicio/terminación."
  },
  {
    id: "codigo_genetico",
    nombre: "Código genético",
    bloque: "bloque5",
    descripcion:
      "Correspondencia entre codones y aminoácidos (y señales start/stop) usada para traducir ARNm a proteína."
  },
  {
    id: "degeneracion_codigo",
    nombre: "Degeneración del código genético",
    bloque: "bloque5",
    descripcion:
      "Propiedad por la cual múltiples codones distintos pueden codificar el mismo aminoácido."
  },
  {
    id: "codon_inicio",
    nombre: "Codón de inicio (start)",
    bloque: "bloque5",
    descripcion:
      "Codón que inicia la traducción; típicamente AUG (Metionina) en el código estándar."
  },
  {
    id: "codon_terminacion",
    nombre: "Codón de terminación (stop)",
    bloque: "bloque5",
    descripcion:
      "Codones que señalan el fin de la traducción; típicamente UAA, UAG o UGA."
  },
  {
    id: "arn_transferencia",
    nombre: "ARN de transferencia (ARNt)",
    bloque: "bloque5",
    descripcion:
      "ARN adaptador que transporta un aminoácido específico y reconoce un codón del ARNm mediante su anticodón."
  },
  {
    id: "anticodon",
    nombre: "Anticodón",
    bloque: "bloque5",
    descripcion:
      "Triplete en ARNt complementario al codón del ARNm; asegura incorporación del aminoácido correcto."
  },
  {
    id: "dogma_central",
    nombre: "Dogma central",
    bloque: "bloque5",
    descripcion:
      "Principio que describe el flujo principal de información genética: ADN → ARN → proteína."
  },

  // === Bloque 6: Tecnologías ADN Recombinante (14) ===
  {
    id: "pcr",
    nombre: "PCR",
    bloque: "bloque6",
    descripcion:
      "Reacción en Cadena de la Polimerasa: técnica de amplificación de un fragmento específico de ADN mediante ciclos de desnaturalización, alineamiento de cebadores y extensión."
  },
  {
    id: "cebador",
    nombre: "Cebador/Primer",
    bloque: "bloque6",
    descripcion:
      "Oligonucleótido corto que se aparea con el molde y provee un extremo 3’ para que la polimerasa inicie síntesis."
  },
  {
    id: "desnaturalizacion",
    nombre: "Desnaturalización (PCR)",
    bloque: "bloque6",
    descripcion:
      "Separación de hebras de ADN por calentamiento para obtener moldes monocatenarios."
  },
  {
    id: "alineamiento_priming",
    nombre: "Alineamiento/Priming (PCR)",
    bloque: "bloque6",
    descripcion:
      "Unión de cebadores a secuencias complementarias flanqueantes al enfriar la reacción."
  },
  {
    id: "extension_pcr",
    nombre: "Extensión (PCR)",
    bloque: "bloque6",
    descripcion:
      "Síntesis de nuevas hebras complementarias a partir de cebadores mediante ADN polimerasa."
  },
  {
    id: "clonacion_molecular",
    nombre: "Clonación molecular",
    bloque: "bloque6",
    descripcion:
      "Inserción de un fragmento de ADN en un vector replicable para producir muchas copias en un organismo hospedero."
  },
  {
    id: "vector_clonacion",
    nombre: "Vector de clonación",
    bloque: "bloque6",
    descripcion:
      "Molécula de ADN (p. ej., plasmidio/virus) capaz de replicarse y transportar un inserto de ADN."
  },
  {
    id: "biblioteca_clones",
    nombre: "Biblioteca de clones",
    bloque: "bloque6",
    descripcion:
      "Colección de clones que representan fragmentos (frecuentemente aleatorios) de un genoma."
  },
  {
    id: "enzima_restriccion",
    nombre: "Enzima de restricción",
    bloque: "bloque6",
    descripcion:
      "Proteína que corta ADN en secuencias específicas (sitios de reconocimiento), generando fragmentos."
  },
  {
    id: "sitio_reconocimiento",
    nombre: "Sitio de reconocimiento",
    bloque: "bloque6",
    descripcion:
      "Secuencia específica (a menudo palindrómica) donde una enzima de restricción se une y corta."
  },
  {
    id: "extremos_romos",
    nombre: "Extremos romos (blunt ends)",
    bloque: "bloque6",
    descripcion:
      "Resultado de un corte que deja extremos sin salientes monocatenarios."
  },
  {
    id: "extremos_pegajosos",
    nombre: "Extremos pegajosos (sticky ends)",
    bloque: "bloque6",
    descripcion:
      "Extremos con salientes monocatenarios complementarios que facilitan unión por hibridación."
  },
  {
    id: "hibridacion",
    nombre: "Hibridación",
    bloque: "bloque6",
    descripcion:
      "Unión de hebras complementarias de ácidos nucleicos por apareamiento de bases."
  },
  {
    id: "ligacion",
    nombre: "Ligación",
    bloque: "bloque6",
    descripcion:
      "Formación de enlaces covalentes para “sellar” el esqueleto azúcar-fosfato y unir fragmentos de ADN."
  },

  // === Bloque 7: Técnicas de Análisis (3) ===
  {
    id: "electroforesis_gel",
    nombre: "Electroforesis en gel",
    bloque: "bloque7",
    descripcion:
      "Técnica que separa fragmentos de ADN por tamaño al migrar en un gel bajo campo eléctrico."
  },
  {
    id: "sonda",
    nombre: "Sonda (probe)",
    bloque: "bloque7",
    descripcion:
      "Oligonucleótido marcado de secuencia conocida que se une por hibridación para detectar una secuencia complementaria."
  },
  {
    id: "microarreglo",
    nombre: "Microarreglo / DNA array",
    bloque: "bloque7",
    descripcion:
      "Superficie con muchas sondas inmovilizadas usada para detectar presencia o abundancia de transcritos (expresión génica) por hibridación."
  },

  // === Bloque 8: Variación y Evolución (7) ===
  {
    id: "variacion_genetica",
    nombre: "Variación genética intraespecífica",
    bloque: "bloque8",
    descripcion:
      "Diferencias de secuencia entre individuos de una misma especie que sustentan diversidad de rasgos."
  },
  {
    id: "genoma_referencia",
    nombre: "Genoma de referencia",
    bloque: "bloque8",
    descripcion:
      "Secuencia representativa (“maestra”) usada como patrón para una especie, aunque individuos difieran en posiciones específicas."
  },
  {
    id: "conservacion_genetica",
    nombre: "Conservación genética",
    bloque: "bloque8",
    descripcion:
      "Presencia de genes o secuencias similares entre especies, indicativa de origen común o función preservada."
  },
  {
    id: "evolucion",
    nombre: "Evolución",
    bloque: "bloque8",
    descripcion:
      "Proceso de cambio heredable en poblaciones a través del tiempo, reflejado en variaciones de secuencia genómica."
  },
  {
    id: "seleccion_natural",
    nombre: "Selección natural",
    bloque: "bloque8",
    descripcion:
      "Mecanismo por el cual variantes genéticas que mejoran el éxito reproductivo tienden a aumentar su frecuencia."
  },
  {
    id: "adaptacion",
    nombre: "Adaptación",
    bloque: "bloque8",
    descripcion:
      "Incremento de frecuencia de rasgos heredables que mejoran el ajuste de una población a su ambiente."
  },
  {
    id: "especiacion",
    nombre: "Especiación",
    bloque: "bloque8",
    descripcion:
      "Formación de nuevas especies cuando poblaciones divergen hasta volverse reproductivamente incompatibles."
  },

  // === Bloque 9: Bioinformática (5) ===
  {
    id: "alfabeto_molecular",
    nombre: "Alfabeto molecular",
    bloque: "bloque9",
    descripcion:
      "Representación secuencial de macromoléculas como “cadenas” (ADN/ARN: 4 letras; proteínas: 20 letras), base conceptual para análisis computacional de secuencias."
  },
  {
    id: "genomica_comparativa",
    nombre: "Genómica comparativa",
    bloque: "bloque9",
    descripcion:
      "Análisis comparado de genomas para identificar similitudes/diferencias, inferir funciones y relaciones evolutivas."
  },
  {
    id: "alineamiento_secuencias",
    nombre: "Alineamiento de secuencias",
    bloque: "bloque9",
    descripcion:
      "Método computacional para comparar secuencias y localizar regiones homólogas o conservadas."
  },
  {
    id: "blast",
    nombre: "BLAST",
    bloque: "bloque9",
    descripcion:
      "Familia de algoritmos rápidos para buscar similitudes entre secuencias y evaluar su significancia estadística."
  },
  {
    id: "bioinformatica",
    nombre: "Bioinformática",
    bloque: "bloque9",
    descripcion:
      "Disciplina que utiliza algoritmos, estadística y modelos computacionales para analizar datos biológicos y decodificar patrones funcionales y evolutivos."
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


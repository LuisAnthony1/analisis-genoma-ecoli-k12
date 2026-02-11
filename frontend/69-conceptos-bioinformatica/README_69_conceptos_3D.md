# Guia / Prompt: Transformar 69 Conceptos de Bioinformatica a una Experiencia 3D Interactiva

## Estructura Actual del Archivo `69_conceptos_bioinformatica.html`

### Datos (todo esta en un solo archivo HTML)

1. **`conceptos`** - Objeto JS con 9 bloques tematicos:
   - `bloque1`: Fundamentos Celulares (4 conceptos, color emerald)
   - `bloque2`: Macromoleculas Informacionales (9 conceptos, color cyan)
   - `bloque3`: Organizacion Genomica (7 conceptos, color purple)
   - `bloque4`: Celulas y Organizacion (4 conceptos, color amber)
   - `bloque5`: Flujo de Informacion Genetica (16 conceptos, color red)
   - `bloque6`: Tecnologias ADN Recombinante (14 conceptos, color pink)
   - `bloque7`: Tecnicas de Analisis (3 conceptos, color blue)
   - `bloque8`: Variacion y Evolucion (7 conceptos, color slate)
   - `bloque9`: Bioinformatica (5 conceptos, color orange)

2. **Cada concepto tiene**: `id` (1-69), `nombre`, y array de `conexiones` (IDs de conceptos relacionados)

3. **`definiciones`** - Objeto con clave=id y valor=texto largo explicando cada concepto

4. **`explicacionesConexion`** - Objeto con clave="id1-id2" (menor-mayor) y valor=texto explicando la relacion entre esos dos conceptos

### UI Actual
- Pagina 2D con cards por bloque, chips por concepto
- Modal al hacer click que muestra definicion + conexiones
- Tema claro/oscuro con switch
- Navegacion anterior/siguiente entre conceptos
- Usa TailwindCSS via CDN, fuentes Google (JetBrains Mono, Space Grotesk)

---

## PROMPT / Instrucciones para Crear la Version 3D

### Objetivo
Crear una visualizacion 3D tipo grafo/red donde los 69 conceptos sean nodos flotando en un espacio 3D, con lineas/aristas animadas conectando conceptos relacionados. El usuario puede rotar, hacer zoom, y hacer click en nodos para explorar.

### Tecnologia Recomendada
- **Three.js** para el renderizado 3D (via CDN: `https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js`)
- **CSS3** para UI superpuesta (paneles de info, controles)
- **Vanilla JS** para logica (mantener sin frameworks)
- Alternativa: **Force-Graph-3D** (`https://unpkg.com/3d-force-graph`) que ya trae fisicas de grafo incluidas

### Estructura del Grafo 3D

```
Nodos (esferas/puntos 3D):
- Cada concepto = 1 nodo esfera
- Color del nodo = color del bloque al que pertenece
- Tamano del nodo = proporcional al numero de conexiones que tiene
- Texto flotante encima del nodo con el nombre del concepto

Aristas (lineas 3D):
- Cada conexion entre conceptos = 1 linea entre nodos
- Color de la linea = gradiente entre los colores de los 2 bloques conectados
- Animacion: particulas/puntos de luz viajando por la linea (como pulsos de informacion)

Clusters:
- Los nodos del mismo bloque deben agruparse naturalmente cerca
- Usar simulacion de fuerzas (force-directed layout) para posicionamiento automatico
- Los bloques deben estar separados pero los conceptos con muchas conexiones inter-bloque deben estar en el borde de su cluster
```

### Interacciones Requeridas

```
1. HOVER sobre nodo:
   - El nodo brilla/crece ligeramente
   - Se iluminan SOLO las aristas conectadas a ese nodo
   - Las demas aristas se atenuan (opacity baja)
   - Tooltip flotante con: nombre + bloque + numero de conexiones

2. CLICK en nodo:
   - Panel lateral derecho se abre con:
     * Nombre del concepto (titulo)
     * Bloque al que pertenece (badge con color)
     * Definicion completa (del objeto `definiciones`)
     * Lista de conexiones con explicacion (del objeto `explicacionesConexion`)
     * Cada conexion es clickeable -> navega al otro nodo
   - La camara hace zoom/orbita suavemente hacia el nodo seleccionado
   - Los nodos conectados pulsan/brillan

3. CLICK en conexion (desde el panel):
   - La camara viaja animadamente del nodo actual al nodo destino
   - Efecto visual de "viaje" por la arista (particula recorriendo la linea)
   - Se abre el panel del nuevo concepto

4. Controles de camara:
   - Click + drag = rotar
   - Scroll = zoom
   - Click derecho + drag = pan
   - Boton "reset camara" para volver a vista general

5. Filtros por bloque:
   - Barra lateral izquierda con los 9 bloques como toggles
   - Cada toggle muestra/oculta los nodos de ese bloque
   - Color coding consistente con el archivo actual
```

### Animaciones Especificas

```
1. CARGA INICIAL:
   - Los nodos aparecen uno a uno desde el centro con efecto de explosion suave
   - Las aristas se dibujan progresivamente despues de que los nodos estan en posicion
   - Duracion total: ~3 segundos

2. ARISTAS ACTIVAS:
   - Particulas de luz (pequenas esferas brillantes) viajan por las aristas continuamente
   - Velocidad lenta, efecto sutil tipo "flujo de datos"
   - Color de la particula = color del nodo origen

3. IDLE (sin interaccion):
   - Rotacion automatica lenta del grafo completo
   - Los nodos flotan/oscilan ligeramente (efecto breathing)
   - Las particulas en aristas siguen fluyendo

4. TRANSICION entre conceptos:
   - Cuando el usuario navega de concepto A a concepto B:
     * Una particula grande y brillante viaja por la arista A->B
     * La camara sigue esa particula suavemente
     * Al llegar a B, efecto de "onda expansiva" sutil
   - Duracion: ~1.5 segundos

5. FONDO:
   - Espacio oscuro con particulas sutiles tipo estrellas
   - Grid sutil en el fondo (similar al actual)
   - Efecto de profundidad/niebla para nodos lejanos
```

### Layout del HTML

```
+--------------------------------------------------+
| NAV: titulo + toggle tema + stats                |
+--------+---------------------------------+-------+
|        |                                 |       |
| FILTROS|      CANVAS 3D (Three.js)       | PANEL |
| POR    |                                 | INFO  |
| BLOQUE |    [grafo 3D interactivo]       | DEL   |
| (izq)  |                                 |CONCEPTO|
|        |                                 | (der) |
|        |                                 |       |
+--------+---------------------------------+-------+
| FOOTER: controles camara + creditos              |
+--------------------------------------------------+
```

### Reutilizacion de Datos

Copiar tal cual del archivo actual estos objetos JS:
- `conceptos` (estructura de bloques, IDs, conexiones)
- `definiciones` (textos descriptivos)
- `explicacionesConexion` (textos de relacion entre conceptos)

NO hay que cambiar los datos, solo la visualizacion.

### Efectos CSS Adicionales

```css
/* Panel de info con efecto glassmorphism */
.panel-info {
    background: rgba(0, 0, 0, 0.7);
    backdrop-filter: blur(20px);
    border: 1px solid rgba(16, 185, 129, 0.3);
    border-radius: 16px;
}

/* Nodo seleccionado - glow effect */
.nodo-activo {
    box-shadow: 0 0 30px currentColor, 0 0 60px currentColor;
}

/* Transicion suave del panel */
.panel-slide {
    transition: transform 0.4s cubic-bezier(0.4, 0, 0.2, 1);
}

/* Texto de concepto flotante sobre nodo 3D - usar CSS2DRenderer de Three.js */
.label-concepto {
    font-family: 'JetBrains Mono', monospace;
    font-size: 11px;
    padding: 2px 8px;
    background: rgba(0, 0, 0, 0.6);
    border-radius: 4px;
    color: white;
    pointer-events: none;
    white-space: nowrap;
}
```

### Ejemplo de Estructura de Codigo Three.js

```javascript
// Pseudocodigo - estructura general
// 1. Crear escena, camara, renderer
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, w/h, 0.1, 1000);
const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });

// 2. Crear nodos como esferas
Object.values(conceptos).forEach(bloque => {
    bloque.conceptos.forEach(concepto => {
        const geometry = new THREE.SphereGeometry(radius, 32, 32);
        const material = new THREE.MeshPhongMaterial({ color: bloque.color, emissive: bloque.color, emissiveIntensity: 0.3 });
        const sphere = new THREE.Mesh(geometry, material);
        // Posicionar con force-directed layout
        scene.add(sphere);
    });
});

// 3. Crear aristas como lineas
conceptos.forEach(concepto => {
    concepto.conexiones.forEach(targetId => {
        const geometry = new THREE.BufferGeometry().setFromPoints([posA, posB]);
        const line = new THREE.Line(geometry, material);
        scene.add(line);
    });
});

// 4. Raycaster para deteccion de clicks/hover
const raycaster = new THREE.Raycaster();
// En cada frame: detectar intersecciones con nodos

// 5. OrbitControls para camara
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.autoRotate = true;
controls.autoRotateSpeed = 0.5;

// 6. Loop de animacion
function animate() {
    requestAnimationFrame(animate);
    // Actualizar particulas en aristas
    // Actualizar fisicas del grafo
    // Animar nodos flotantes
    controls.update();
    renderer.render(scene, camera);
}
```

### Alternativa Mas Simple: 3d-force-graph

Si Three.js puro es muy complejo, usar la libreria `3d-force-graph` que maneja automaticamente:
- Posicionamiento por fuerzas
- Renderizado 3D
- Controles de camara
- Hover/click en nodos

```javascript
import ForceGraph3D from '3d-force-graph';

const graph = ForceGraph3D()(document.getElementById('graph-container'))
    .graphData({
        nodes: todosLosConceptos.map(c => ({ id: c.id, name: c.nombre, group: c.bloqueId, color: c.color })),
        links: todasLasConexiones.map(([source, target]) => ({ source, target }))
    })
    .nodeColor(node => node.color)
    .nodeLabel(node => node.name)
    .onNodeClick(node => mostrarPanel(node));
```

### Requisitos Finales

- Todo en UN solo archivo HTML (como el actual)
- Librerias via CDN (no npm/build)
- Responsive (que funcione en movil con controles tactiles)
- Performance: debe correr fluido con 69 nodos y ~200 aristas
- Mantener el toggle de tema claro/oscuro
- Mantener los mismos colores por bloque
- Los datos (conceptos, definiciones, explicaciones) son IDENTICOS al archivo actual

---

**Autor**: Luis Anthony Mamani Mescco
**Curso**: Bioinformatica 2025-2V // UNSAAC
**Archivo fuente**: `69_conceptos_bioinformatica.html`

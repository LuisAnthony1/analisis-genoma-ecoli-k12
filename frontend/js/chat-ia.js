/**
 * GenomeHub - Chat IA Especialista
 *
 * Modal de chat conectado a Gemini API via servidor.py
 * Permite conversar sobre los resultados de analisis genomico.
 */

const ChatIA = {
    modal: null,
    historial: [],
    analysisType: '',
    genome: '',
    contextData: null,

    open(analysisType, genome) {
        if (!genome) {
            showNotification('Selecciona un genoma primero', 'warning');
            return;
        }

        this.analysisType = analysisType || 'genes';

        // Si es el mismo genoma y el modal ya existe, solo mostrarlo
        if (this.genome === genome && this.modal) {
            this.modal.style.display = 'flex';
            setTimeout(() => document.getElementById('chat-ia-input')?.focus(), 100);
            return;
        }

        // Nuevo genoma o primera vez: resetear todo
        this.genome = genome;
        this.historial = [];
        this.contextData = null;

        this._crearModal();
        this._cargarContexto();
    },

    _crearModal() {
        // Remover modal anterior si existe
        if (this.modal) this.modal.remove();

        const modal = document.createElement('div');
        modal.id = 'chat-ia-modal';
        modal.className = 'fixed inset-0 z-50 flex items-center justify-center p-4 bg-black/50 backdrop-blur-sm';
        modal.onclick = (e) => { if (e.target === modal) this.close(); };

        modal.innerHTML = `
            <div class="bg-card rounded-2xl w-full max-w-2xl max-h-[85vh] flex flex-col border border-slate-200 shadow-2xl" onclick="event.stopPropagation()">
                <!-- Header -->
                <div class="flex items-center justify-between px-6 py-4 border-b border-slate-200 dark:border-slate-700">
                    <div class="flex items-center gap-3">
                        <div class="w-10 h-10 rounded-xl bg-gradient-to-br from-violet-500 to-fuchsia-500 flex items-center justify-center text-xl shadow-lg">
                            ✨
                        </div>
                        <div>
                            <h3 class="font-bold text-primary">Especialista en Genomica</h3>
                            <p class="text-xs text-secondary">Analizando: ${this.genome.replace(/_/g, ' ')} - ${this.analysisType}</p>
                        </div>
                    </div>
                    <button onclick="ChatIA.close()" class="w-8 h-8 flex items-center justify-center rounded-lg hover:bg-slate-100 dark:hover:bg-slate-800 text-secondary transition">
                        ✕
                    </button>
                </div>

                <!-- Mensajes -->
                <div id="chat-ia-messages" class="flex-1 overflow-y-auto px-6 py-4 space-y-4" style="min-height: 300px; max-height: 50vh;">
                    <div class="flex gap-3">
                        <div class="w-8 h-8 rounded-lg bg-gradient-to-br from-violet-500 to-fuchsia-500 flex items-center justify-center text-sm flex-shrink-0">✨</div>
                        <div class="bg-slate-50 dark:bg-slate-800 rounded-xl rounded-tl-sm px-4 py-3 max-w-[85%]">
                            <p class="text-sm text-primary">Hola! Soy tu especialista en genomica. Estoy analizando los resultados de <strong>${this.analysisType}</strong> del genoma <strong>${this.genome.replace(/_/g, ' ')}</strong>.</p>
                            <p class="text-sm text-primary mt-2">Puedes preguntarme sobre:</p>
                            <ul class="text-sm text-secondary mt-1 space-y-1 list-disc list-inside">
                                <li>Interpretacion de los resultados</li>
                                <li>Comparacion con otros organismos</li>
                                <li>Significado biologico de las metricas</li>
                                <li>Cualquier duda sobre genomica</li>
                            </ul>
                        </div>
                    </div>
                    <div id="chat-ia-loading-context" class="flex gap-3">
                        <div class="w-8 h-8 rounded-lg bg-slate-200 dark:bg-slate-700 flex items-center justify-center text-sm flex-shrink-0">⏳</div>
                        <div class="bg-slate-50 dark:bg-slate-800 rounded-xl px-4 py-3">
                            <p class="text-xs text-secondary">Cargando datos del analisis...</p>
                        </div>
                    </div>
                </div>

                <!-- Quick Actions -->
                <div class="px-6 py-2 border-t border-slate-100 dark:border-slate-800 flex gap-2 flex-wrap">
                    <button onclick="ChatIA._quickAction('crispr')" class="px-3 py-1.5 bg-violet-500/10 text-violet-500 text-xs font-medium rounded-lg hover:bg-violet-500/20 transition">
                        Sitios CRISPR-Cas9
                    </button>
                    <button onclick="ChatIA._quickAction('mutacion')" class="px-3 py-1.5 bg-amber-500/10 text-amber-500 text-xs font-medium rounded-lg hover:bg-amber-500/20 transition">
                        Mutaciones clave
                    </button>
                    <button onclick="ChatIA._quickAction('resumen')" class="px-3 py-1.5 bg-emerald-500/10 text-emerald-500 text-xs font-medium rounded-lg hover:bg-emerald-500/20 transition">
                        Resumen completo
                    </button>
                    <button onclick="ChatIA._quickAction('patogenicidad')" class="px-3 py-1.5 bg-red-500/10 text-red-500 text-xs font-medium rounded-lg hover:bg-red-500/20 transition">
                        Patogenicidad
                    </button>
                </div>

                <!-- Input -->
                <div class="px-6 py-4 border-t border-slate-200 dark:border-slate-700">
                    <div class="flex gap-3">
                        <input
                            type="text"
                            id="chat-ia-input"
                            placeholder="Escribe tu pregunta..."
                            class="flex-1 px-4 py-3 bg-slate-50 dark:bg-slate-800 rounded-xl border border-slate-200 dark:border-slate-700 focus:border-violet-500 focus:ring-2 focus:ring-violet-500/20 outline-none text-sm text-primary"
                            onkeydown="if(event.key==='Enter') ChatIA.enviar()"
                        >
                        <button
                            id="chat-ia-send-btn"
                            onclick="ChatIA.enviar()"
                            class="px-5 py-3 bg-gradient-to-r from-violet-500 to-fuchsia-500 hover:from-violet-600 hover:to-fuchsia-600 text-white rounded-xl font-medium text-sm transition shadow-lg"
                        >
                            Enviar
                        </button>
                    </div>
                    <p class="text-xs text-secondary mt-2 text-center">Potenciado por Gemini IA</p>
                </div>
            </div>
        `;

        document.body.appendChild(modal);
        this.modal = modal;

        // Focus en input
        setTimeout(() => document.getElementById('chat-ia-input')?.focus(), 100);
    },

    async _cargarContexto() {
        // Cargar TODOS los archivos JSON de resultados para dar contexto completo a la IA
        const archivosBase = [
            `analisis_genes_${this.genome}.json`,
            `analisis_codones_${this.genome}.json`,
            `analisis_distancias_${this.genome}.json`
        ];

        // Buscar comparaciones disponibles
        try {
            const resp = await fetch(`/api/results?type=tablas&genome=${this.genome}`);
            const listado = await resp.json();
            if (listado.success && listado.results) {
                for (const f of listado.results) {
                    if (f.extension === 'json' && !archivosBase.includes(f.filename)) {
                        archivosBase.push(f.filename);
                    }
                }
            }
        } catch (e) {
            console.error('Error listando archivos:', e);
        }

        const resumen = {};
        let algoCargado = false;

        // Cargar todos los JSON en paralelo
        const promesas = archivosBase.map(async (filename) => {
            try {
                const data = await DashboardRenderer.fetchResultData(this.genome, filename);
                if (data) {
                    const key = filename.replace(`_${this.genome}`, '').replace('.json', '');
                    resumen[key] = this._resumirDatos(data, filename);
                    algoCargado = true;
                }
            } catch (e) {
                // Archivo no existe, ignorar
            }
        });

        await Promise.all(promesas);

        if (algoCargado) {
            this.contextData = resumen;
        }

        this._removerLoadingContext();
    },

    _resumirDatos(data, filename = '') {
        // Crear resumen conciso segun tipo de archivo
        const resumen = {};

        // Genes
        if (data.estadisticas_generales) {
            resumen.estadisticas_generales = data.estadisticas_generales;
        }
        if (data.comparacion_literatura && Object.keys(data.comparacion_literatura).length > 0) {
            resumen.comparacion_literatura = data.comparacion_literatura;
        }
        if (data.distribucion_tamanos) resumen.distribucion_tamanos = data.distribucion_tamanos;
        if (data.genes_extremos) resumen.genes_extremos = data.genes_extremos;
        if (data.analisis_genes_vs_cds) resumen.genes_vs_cds = data.analisis_genes_vs_cds;

        // Codones
        if (data.conteo_64_codones) {
            // Solo top 10 codones para no exceder tokens
            const top = data.conteo_64_codones.tabla_codones;
            if (top) resumen.top_codones = top.slice(0, 10);
        }
        if (data.codones_inicio) resumen.codones_inicio = data.codones_inicio;
        if (data.codones_parada) resumen.codones_parada = data.codones_parada;
        if (data.contenido_gc) resumen.contenido_gc = data.contenido_gc;

        // Distancias
        if (data.estadisticas_generales && data.distribucion_tipos) {
            resumen.distribucion_distancias = data.distribucion_tipos;
        }

        // Comparacion de genomas
        if (data.organismos_comparados) resumen.organismos_comparados = data.organismos_comparados;
        if (data.metricas_generales) resumen.metricas_generales = data.metricas_generales;
        if (data.genes_virulencia) {
            // Resumir virulencia (solo totales, no lista completa)
            const vir = {};
            for (const [k, v] of Object.entries(data.genes_virulencia)) {
                if (typeof v === 'object' && v.total !== undefined) {
                    vir[k] = { total: v.total, categorias: v.categorias || v.por_categoria };
                } else {
                    vir[k] = v;
                }
            }
            resumen.genes_virulencia = vir;
        }
        if (data.resumen_interpretativo) resumen.resumen = data.resumen_interpretativo;
        if (data.uso_codones) resumen.uso_codones = data.uso_codones;
        if (data.distribucion_tamanos) resumen.distribucion_tamanos = data.distribucion_tamanos;
        if (data.distribucion_gc) resumen.distribucion_gc = data.distribucion_gc;
        if (data.interpretacion_ia) resumen.interpretacion_ia = data.interpretacion_ia;

        // Si no se extrajo nada relevante, copiar campos de primer nivel (sin arrays largos)
        if (Object.keys(resumen).length === 0) {
            for (const [k, v] of Object.entries(data)) {
                if (Array.isArray(v) && v.length > 20) continue; // Omitir arrays grandes
                resumen[k] = v;
            }
        }

        return resumen;
    },

    _removerLoadingContext() {
        const el = document.getElementById('chat-ia-loading-context');
        if (el) {
            if (this.contextData) {
                el.querySelector('p').textContent = 'Datos cargados. Listo para responder!';
                el.querySelector('div:first-child').textContent = '✅';
            } else {
                el.querySelector('p').textContent = 'No se encontraron datos de analisis.';
                el.querySelector('div:first-child').textContent = '⚠️';
            }
        }
    },

    async enviar() {
        const input = document.getElementById('chat-ia-input');
        const mensaje = input.value.trim();
        if (!mensaje) return;

        input.value = '';
        this._agregarMensaje(mensaje, 'user');
        this._mostrarTyping();

        const btn = document.getElementById('chat-ia-send-btn');
        btn.disabled = true;
        input.disabled = true;

        try {
            const genomeName = this.genome.replace(/_/g, ' ');
            const contexto = this.contextData
                ? `RESULTADOS REALES DEL ANALISIS del genoma ${genomeName}:\n${JSON.stringify(this.contextData, null, 1)}`
                : `Genoma: ${genomeName}, Analisis: ${this.analysisType} (sin datos cargados)`;

            const systemPrompt = `Eres un especialista en genomica y bioinformatica. Responde SIEMPRE en espanol.
Estas analizando resultados REALES de un analisis genomico de "${genomeName}".

IMPORTANTE: Se te proporcionan datos REALES del analisis. USA ESOS DATOS para responder con valores especificos y exactos.
- Cuando te pregunten "cuantos genes?", responde con el numero exacto de los datos (ej: total_genes)
- Cuando te pregunten sobre GC%, da el valor exacto de los resultados
- Interpreta los resultados biologicamente: que significan, si son normales, que implicaciones tienen
- Si hay datos de comparacion, usalos para contextualizar
- Si algo es patogeno o tiene genes de virulencia, explica las implicaciones
- Se conciso pero informativo. Maximo 3-4 parrafos por respuesta.
- No repitas "basandome en los datos proporcionados" cada vez, simplemente responde con los datos.`;

            const resp = await fetch('/api/chat', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    message: mensaje,
                    context: contexto,
                    system: systemPrompt,
                    genome: this.genome
                })
            });

            const data = await resp.json();

            this._removerTyping();

            if (data.success && data.response) {
                this._agregarMensaje(data.response, 'ia');
            } else {
                this._agregarMensaje(data.error || 'Error al obtener respuesta', 'error');
            }
        } catch (error) {
            this._removerTyping();
            this._agregarMensaje('Error de conexion: ' + error.message, 'error');
        } finally {
            btn.disabled = false;
            input.disabled = false;
            input.focus();
        }
    },

    _agregarMensaje(texto, tipo) {
        const container = document.getElementById('chat-ia-messages');
        const div = document.createElement('div');
        div.className = 'flex gap-3' + (tipo === 'user' ? ' flex-row-reverse' : '');

        const avatar = tipo === 'user' ? '👤' : tipo === 'error' ? '⚠️' : '✨';
        const avatarBg = tipo === 'user'
            ? 'bg-emerald-500'
            : tipo === 'error' ? 'bg-red-500' : 'bg-gradient-to-br from-violet-500 to-fuchsia-500';
        const bubbleBg = tipo === 'user'
            ? 'bg-emerald-500 text-white rounded-tr-sm'
            : tipo === 'error' ? 'bg-red-500/10 text-red-500 rounded-tl-sm' : 'bg-slate-50 dark:bg-slate-800 text-primary rounded-tl-sm';

        // Formatear texto IA con markdown basico
        let htmlTexto = texto;
        if (tipo === 'ia') {
            htmlTexto = this._formatearMarkdown(texto);
        }

        div.innerHTML = `
            <div class="w-8 h-8 rounded-lg ${avatarBg} flex items-center justify-center text-sm flex-shrink-0 text-white">${avatar}</div>
            <div class="${bubbleBg} rounded-xl px-4 py-3 max-w-[85%]">
                <div class="text-sm leading-relaxed">${htmlTexto}</div>
            </div>
        `;

        container.appendChild(div);
        container.scrollTop = container.scrollHeight;

        this.historial.push({ tipo, texto });
    },

    _formatearMarkdown(texto) {
        return texto
            .replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>')
            .replace(/\*(.*?)\*/g, '<em>$1</em>')
            .replace(/`(.*?)`/g, '<code class="bg-slate-200 dark:bg-slate-700 px-1 rounded text-xs">$1</code>')
            .replace(/\n- /g, '\n• ')
            .replace(/\n/g, '<br>');
    },

    _mostrarTyping() {
        const container = document.getElementById('chat-ia-messages');
        const div = document.createElement('div');
        div.id = 'chat-ia-typing';
        div.className = 'flex gap-3';
        div.innerHTML = `
            <div class="w-8 h-8 rounded-lg bg-gradient-to-br from-violet-500 to-fuchsia-500 flex items-center justify-center text-sm flex-shrink-0">✨</div>
            <div class="bg-slate-50 dark:bg-slate-800 rounded-xl px-4 py-3">
                <div class="flex gap-1">
                    <div class="w-2 h-2 bg-violet-400 rounded-full animate-bounce" style="animation-delay: 0ms"></div>
                    <div class="w-2 h-2 bg-violet-400 rounded-full animate-bounce" style="animation-delay: 150ms"></div>
                    <div class="w-2 h-2 bg-violet-400 rounded-full animate-bounce" style="animation-delay: 300ms"></div>
                </div>
            </div>
        `;
        container.appendChild(div);
        container.scrollTop = container.scrollHeight;
    },

    _removerTyping() {
        const el = document.getElementById('chat-ia-typing');
        if (el) el.remove();
    },

    _quickAction(type) {
        const genomeName = this.genome.replace(/_/g, ' ');
        const prompts = {
            crispr: `Basandote en los datos del genoma ${genomeName}, sugiere 3-5 posibles sitios CRISPR-Cas9 para modificacion genomica. Para cada sitio indica:\n1. Gen objetivo y su funcion\n2. Secuencia guia (sgRNA) sugerida de 20nt + PAM (NGG)\n3. Justificacion biologica\n4. Efectos fenotipicos esperados de la mutacion\nConsidera genes de virulencia, metabolismo y resistencia.`,
            mutacion: `Analiza los genes mas importantes del genoma ${genomeName} y describe las implicaciones biologicas de mutaciones en los 5 genes mas criticos. Para cada gen: funcion normal, que pasaria si se inactiva (knockout), y aplicaciones potenciales en biotecnologia.`,
            resumen: `Genera un resumen ejecutivo completo del analisis genomico de ${genomeName}. Incluye: hallazgos principales, metricas clave, comparacion con organismos similares, y conclusiones biologicas relevantes. Usa los datos reales proporcionados.`,
            patogenicidad: `Analiza la patogenicidad de ${genomeName} basandote en los datos del analisis. Evalua: genes de virulencia encontrados, factores de patogenicidad, islas de patogenicidad potenciales, y riesgo biologico. Compara con organismos patogenos conocidos si hay datos de comparacion.`
        };

        const input = document.getElementById('chat-ia-input');
        if (input && prompts[type]) {
            input.value = prompts[type];
            this.enviar();
        }
    },

    close() {
        if (this.modal) {
            this.modal.style.display = 'none';
        }
    },

    destroy() {
        if (this.modal) {
            this.modal.remove();
            this.modal = null;
        }
        this.historial = [];
        this.contextData = null;
        this.genome = '';
    }
};

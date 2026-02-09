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
                    <p class="text-xs text-secondary mt-2 text-center">Powered by Gemini AI</p>
                </div>
            </div>
        `;

        document.body.appendChild(modal);
        this.modal = modal;

        // Focus en input
        setTimeout(() => document.getElementById('chat-ia-input')?.focus(), 100);
    },

    async _cargarContexto() {
        const fileMap = {
            genes: `analisis_genes_${this.genome}.json`,
            codones: `analisis_codones_${this.genome}.json`,
            distancias: `analisis_distancias_${this.genome}.json`
        };

        const filename = fileMap[this.analysisType];
        if (!filename) {
            this._removerLoadingContext();
            return;
        }

        try {
            const data = await DashboardRenderer.fetchResultData(this.genome, filename);
            if (data) {
                // Resumir datos para no exceder limite de tokens
                this.contextData = this._resumirDatos(data);
            }
        } catch (e) {
            console.error('Error cargando contexto:', e);
        }

        this._removerLoadingContext();
    },

    _resumirDatos(data) {
        // Crear resumen conciso de los datos para enviar como contexto
        const resumen = {};
        const copiarCampos = (obj, campos) => {
            const r = {};
            for (const c of campos) {
                if (obj[c] !== undefined) r[c] = obj[c];
            }
            return r;
        };

        if (data.estadisticas_generales) {
            resumen.estadisticas = data.estadisticas_generales;
        }
        if (data.comparacion_literatura) {
            resumen.literatura = data.comparacion_literatura;
        }
        if (data.codones_inicio) resumen.codones_inicio = data.codones_inicio;
        if (data.codones_parada) resumen.codones_parada = data.codones_parada;
        if (data.contenido_gc) resumen.contenido_gc = data.contenido_gc;
        if (data.distribucion_tipos) resumen.distribucion_distancias = data.distribucion_tipos;
        if (data.genes_extremos) resumen.genes_extremos = data.genes_extremos;

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
            const contexto = this.contextData
                ? `Datos del analisis de ${this.analysisType} del genoma ${this.genome}:\n${JSON.stringify(this.contextData, null, 1)}`
                : `Genoma: ${this.genome}, Analisis: ${this.analysisType}`;

            const systemPrompt = `Eres un especialista en genomica y bioinformatica. Responde en espanol.
Estas analizando resultados de un analisis de ${this.analysisType} del genoma ${this.genome.replace(/_/g, ' ')}.
Responde de forma clara y concisa. Usa datos especificos de los resultados cuando sea relevante.
Si no tienes los datos exactos, indica que los datos no estan disponibles.`;

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

    close() {
        if (this.modal) {
            this.modal.remove();
            this.modal = null;
        }
    }
};

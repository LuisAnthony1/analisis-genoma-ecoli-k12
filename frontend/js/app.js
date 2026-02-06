/**
 * GenomeHub - Core Application
 *
 * Maneja navegación, estado global, temas y utilidades comunes
 */

// =============================================================================
// ESTADO GLOBAL
// =============================================================================

const AppState = {
    currentSection: 'search',
    selectedGenomes: new Set(),
    libraryGenomes: [],
    theme: localStorage.getItem('theme') || 'dark',
    apiBase: window.location.hostname === 'localhost' ? '' : ''
};

// =============================================================================
// NAVEGACIÓN ENTRE SECCIONES
// =============================================================================

function showSection(sectionName) {
    console.log(`[APP] Navegando a sección: ${sectionName}`);

    // Ocultar todas las secciones
    document.querySelectorAll('.section-panel').forEach(section => {
        section.classList.remove('active');
        section.style.display = 'none';
    });

    // Mostrar sección seleccionada
    const targetSection = document.getElementById(`section-${sectionName}`);
    if (targetSection) {
        targetSection.classList.add('active');
        targetSection.style.display = 'block';
    }

    // Actualizar botones de navegación
    document.querySelectorAll('.nav-btn').forEach(btn => {
        btn.classList.remove('active', 'bg-emerald-500', 'text-white');
        btn.classList.add('text-secondary');
    });

    const activeBtn = document.querySelector(`[data-section="${sectionName}"]`);
    if (activeBtn) {
        activeBtn.classList.add('active', 'bg-emerald-500', 'text-white');
        activeBtn.classList.remove('text-secondary');
    }

    AppState.currentSection = sectionName;

    // Cargar datos según la sección
    if (sectionName === 'library') {
        loadLibrary();
    } else if (sectionName === 'results') {
        loadResults('tablas');
    }
}

// =============================================================================
// SISTEMA DE TEMAS
// =============================================================================

function toggleTheme() {
    AppState.theme = AppState.theme === 'dark' ? 'light' : 'dark';
    localStorage.setItem('theme', AppState.theme);
    applyTheme(AppState.theme);
}

function applyTheme(theme) {
    document.body.classList.remove('light', 'dark');
    document.body.classList.add(theme);

    const themeToggle = document.getElementById('theme-toggle');
    const themeSwitch = document.getElementById('themeSwitch');

    if (themeToggle) {
        themeToggle.textContent = theme === 'dark' ? '☀️' : '🌙';
    }

    if (themeSwitch) {
        themeSwitch.className = `theme-switch ${theme}`;
    }
}

// =============================================================================
// NOTIFICACIONES TOAST
// =============================================================================

function showNotification(message, type = 'info') {
    console.log(`[${type.toUpperCase()}] ${message}`);

    const container = document.getElementById('notifications-container') || createNotificationContainer();

    const notification = document.createElement('div');
    notification.className = `notification toast ${type} animate-slide-in`;

    const colors = {
        success: 'bg-emerald-500',
        error: 'bg-red-500',
        warning: 'bg-amber-500',
        info: 'bg-cyan-500'
    };

    const icons = {
        success: '✓',
        error: '✕',
        warning: '!',
        info: 'ⓘ'
    };

    notification.innerHTML = `
        <div class="flex items-center gap-3 px-4 py-3 rounded-lg shadow-lg ${colors[type]} text-white">
            <span class="text-xl font-bold">${icons[type]}</span>
            <span class="flex-1">${message}</span>
            <button onclick="this.parentElement.parentElement.remove()" class="text-white hover:text-gray-200">
                ✕
            </button>
        </div>
    `;

    container.appendChild(notification);

    // Auto-remover después de 5 segundos
    setTimeout(() => {
        notification.style.opacity = '0';
        setTimeout(() => notification.remove(), 300);
    }, 5000);
}

function createNotificationContainer() {
    const container = document.createElement('div');
    container.id = 'notifications-container';
    container.className = 'fixed top-20 right-4 z-50 space-y-2 max-w-md';
    document.body.appendChild(container);
    return container;
}

// =============================================================================
// UTILIDADES
// =============================================================================

function formatFileSize(bytes) {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
}

function formatDate(isoString) {
    if (!isoString) return 'N/A';
    const date = new Date(isoString);
    return date.toLocaleDateString('es-ES', {
        year: 'numeric',
        month: 'short',
        day: 'numeric',
        hour: '2-digit',
        minute: '2-digit'
    });
}

function formatNumber(num) {
    return new Intl.NumberFormat('es-ES').format(num);
}

function showLoading(elementId) {
    const element = document.getElementById(elementId);
    if (element) {
        element.classList.remove('hidden');
    }
}

function hideLoading(elementId) {
    const element = document.getElementById(elementId);
    if (element) {
        element.classList.add('hidden');
    }
}

// =============================================================================
// INICIALIZACIÓN
// =============================================================================

document.addEventListener('DOMContentLoaded', () => {
    console.log('[APP] Inicializando GenomeHub...');

    // Navegación
    document.querySelectorAll('.nav-btn').forEach(btn => {
        btn.addEventListener('click', (e) => {
            const section = e.currentTarget.getAttribute('data-section');
            showSection(section);
        });
    });

    // Toggle tema
    const themeToggle = document.getElementById('theme-toggle');
    if (themeToggle) {
        themeToggle.addEventListener('click', toggleTheme);
    }

    const themeSwitch = document.getElementById('themeSwitch');
    if (themeSwitch) {
        themeSwitch.addEventListener('click', toggleTheme);
    }

    // Aplicar tema inicial
    applyTheme(AppState.theme);

    // Mostrar sección inicial
    showSection('search');

    console.log('[APP] GenomeHub inicializado correctamente');
});

// Estilos CSS para notificaciones (agregar dinámicamente)
const style = document.createElement('style');
style.textContent = `
    .animate-slide-in {
        animation: slideIn 0.3s ease-out;
    }
    @keyframes slideIn {
        from {
            transform: translateX(100%);
            opacity: 0;
        }
        to {
            transform: translateX(0);
            opacity: 1;
        }
    }
    .notification {
        transition: opacity 0.3s ease;
    }
`;
document.head.appendChild(style);

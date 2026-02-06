<?php
/**
 * Configuración global para la API de GenomeHub
 *
 * Este archivo define constantes, funciones de utilidad y configuración
 * para todos los endpoints de la API.
 */

// =============================================================================
// CONFIGURACIÓN DE RUTAS
// =============================================================================

define('BACKEND_DIR', dirname(__DIR__));
define('SCRIPTS_DIR', BACKEND_DIR . DIRECTORY_SEPARATOR . 'scripts');
define('CRUDO_DIR', BACKEND_DIR . DIRECTORY_SEPARATOR . 'crudo');
define('RESULTADOS_DIR', BACKEND_DIR . DIRECTORY_SEPARATOR . 'resultados');

// Detectar entorno (local vs servidor remoto)
$is_local = in_array($_SERVER['SERVER_NAME'] ?? 'localhost', ['localhost', '127.0.0.1', '::1']);
define('IS_LOCAL', $is_local);
define('PYTHON_BIN', $is_local ? 'python' : 'python3');

// =============================================================================
// CONFIGURACIÓN DE PHP
// =============================================================================

// Timeouts largos para descargas y análisis
set_time_limit(300); // 5 minutos
ini_set('max_execution_time', '300');
ini_set('memory_limit', '512M');

// =============================================================================
// HEADERS HTTP
// =============================================================================

// CORS para desarrollo local
header('Access-Control-Allow-Origin: *');
header('Access-Control-Allow-Methods: GET, POST, OPTIONS');
header('Access-Control-Allow-Headers: Content-Type');
header('Content-Type: application/json; charset=utf-8');

// Manejar preflight OPTIONS
if ($_SERVER['REQUEST_METHOD'] === 'OPTIONS') {
    http_response_code(200);
    exit();
}

// =============================================================================
// FUNCIONES DE UTILIDAD
// =============================================================================

/**
 * Sanitiza input del usuario
 */
function sanitize_input($data) {
    if (is_array($data)) {
        return array_map('sanitize_input', $data);
    }
    return htmlspecialchars(strip_tags(trim($data)), ENT_QUOTES, 'UTF-8');
}

/**
 * Valida que un accession ID sea válido (formato NCBI)
 * Ejemplos válidos: NC_000913.3, U00096.3, CP000243.1
 */
function validate_accession_id($id) {
    // Formato: 1-2 letras mayúsculas, opcional guión bajo, 6+ dígitos, punto, versión
    return preg_match('/^[A-Z]{1,2}_?\d{6,}\.\d+$/', $id) === 1;
}

/**
 * Sanitiza nombre de organismo
 * Solo permite: letras, números, espacios, guiones, puntos
 */
function sanitize_organism_name($name) {
    return preg_replace('/[^a-zA-Z0-9\s\-_\.]/', '', $name);
}

/**
 * Valida que una ruta sea segura (no contenga ../)
 */
function is_safe_path($path) {
    $realpath = realpath($path);
    // La ruta debe existir y estar dentro de BACKEND_DIR
    return $realpath !== false && strpos($realpath, BACKEND_DIR) === 0;
}

/**
 * Ejecuta un script Python y retorna el resultado
 *
 * @param string $script Nombre del script (ej: "buscar_ncbi.py")
 * @param array $args Argumentos para el script
 * @param int $timeout Timeout en segundos (default: 30)
 * @return array ['output' => string, 'return_code' => int, 'success' => bool]
 */
function execute_python($script, $args = [], $timeout = 30) {
    $script_path = SCRIPTS_DIR . DIRECTORY_SEPARATOR . $script;

    // Verificar que el script existe
    if (!file_exists($script_path)) {
        return [
            'output' => "Script no encontrado: $script",
            'return_code' => 1,
            'success' => false
        ];
    }

    // Construir comando
    $command = PYTHON_BIN . ' "' . $script_path . '"';

    // Agregar argumentos (escapados)
    foreach ($args as $arg) {
        $command .= ' "' . escapeshellarg($arg) . '"';
    }

    // Redirigir stderr a stdout
    $command .= ' 2>&1';

    // Ejecutar comando
    $output = [];
    $return_code = 0;

    // Nota: En Windows no hay comando 'timeout', usamos exec directo
    // Para servidores Linux, se puede agregar: timeout {$timeout}s antes del comando
    if (!IS_LOCAL && strtoupper(substr(PHP_OS, 0, 3)) !== 'WIN') {
        $command = "timeout {$timeout}s " . $command;
    }

    exec($command, $output, $return_code);

    return [
        'output' => implode("\n", $output),
        'return_code' => $return_code,
        'success' => $return_code === 0
    ];
}

/**
 * Retorna respuesta JSON de error
 */
function json_error($message, $code = 400, $details = null) {
    http_response_code($code);
    $response = [
        'success' => false,
        'error' => $message
    ];
    if ($details !== null) {
        $response['details'] = $details;
    }
    echo json_encode($response, JSON_PRETTY_PRINT | JSON_UNESCAPED_UNICODE);
    exit();
}

/**
 * Retorna respuesta JSON de éxito
 */
function json_success($data = [], $message = null) {
    $response = array_merge(['success' => true], $data);
    if ($message !== null) {
        $response['message'] = $message;
    }
    echo json_encode($response, JSON_PRETTY_PRINT | JSON_UNESCAPED_UNICODE);
    exit();
}

/**
 * Log de debugging (solo en local)
 */
function debug_log($message) {
    if (IS_LOCAL) {
        error_log("[GenomeHub DEBUG] " . print_r($message, true));
    }
}

// =============================================================================
// INICIALIZACIÓN
// =============================================================================

// Crear directorios si no existen
if (!is_dir(CRUDO_DIR)) {
    mkdir(CRUDO_DIR, 0755, true);
}

if (!is_dir(RESULTADOS_DIR . DIRECTORY_SEPARATOR . 'tablas')) {
    mkdir(RESULTADOS_DIR . DIRECTORY_SEPARATOR . 'tablas', 0755, true);
}

if (!is_dir(RESULTADOS_DIR . DIRECTORY_SEPARATOR . 'figuras')) {
    mkdir(RESULTADOS_DIR . DIRECTORY_SEPARATOR . 'figuras', 0755, true);
}

debug_log("Config loaded - Python: " . PYTHON_BIN . " | Local: " . (IS_LOCAL ? 'YES' : 'NO'));

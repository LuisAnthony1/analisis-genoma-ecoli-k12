<?php
/**
 * Endpoint: Obtener resultados de análisis
 *
 * GET /backend/api/get_results.php?type=tablas|figuras
 *
 * Parámetros:
 * - type (string): 'tablas' para JSON/CSV o 'figuras' para PNG (default: tablas)
 *
 * Retorna:
 * - success: true/false
 * - type: tipo de resultados ('tablas' o 'figuras')
 * - count: número de archivos encontrados
 * - results: array de archivos con metadata
 */

require_once 'config.php';

// =============================================================================
// VALIDACIÓN DE PARÁMETROS
// =============================================================================

$type = isset($_GET['type']) ? sanitize_input($_GET['type']) : 'tablas';

if (!in_array($type, ['tablas', 'figuras'])) {
    json_error('Parámetro "type" inválido. Valores válidos: tablas, figuras', 400);
}

debug_log("Getting results: type=$type");

// =============================================================================
// ESCANEAR ARCHIVOS
// =============================================================================

$results = [];

if ($type === 'tablas') {
    // Buscar archivos JSON y CSV
    $tablas_dir = RESULTADOS_DIR . DIRECTORY_SEPARATOR . 'tablas' . DIRECTORY_SEPARATOR;

    $json_files = glob($tablas_dir . '*.json');
    $csv_files = glob($tablas_dir . '*.csv');

    $files = array_merge($json_files ?: [], $csv_files ?: []);

} elseif ($type === 'figuras') {
    // Buscar archivos PNG
    $figuras_dir = RESULTADOS_DIR . DIRECTORY_SEPARATOR . 'figuras' . DIRECTORY_SEPARATOR;

    $files = glob($figuras_dir . '*.png') ?: [];
}

// Verificar errores al escanear
if ($files === false) {
    json_error("Error al escanear directorio de resultados: $type", 500);
}

// =============================================================================
// EXTRAER METADATA DE CADA ARCHIVO
// =============================================================================

foreach ($files as $file) {
    $filename = basename($file);
    $extension = pathinfo($file, PATHINFO_EXTENSION);
    $size_kb = round(filesize($file) / 1024, 2);
    $modified = date('Y-m-d H:i:s', filemtime($file));

    // Ruta relativa desde la raíz del proyecto para acceso web
    $relative_path = 'backend/resultados/' . $type . '/' . $filename;

    $file_info = [
        'filename' => $filename,
        'path' => $relative_path,
        'extension' => $extension,
        'size_kb' => $size_kb,
        'modified' => $modified,
        'timestamp' => filemtime($file)
    ];

    // Para archivos JSON, intentar leer preview de contenido
    if ($extension === 'json' && filesize($file) < 100000) { // Solo si < 100KB
        $content = @file_get_contents($file);
        if ($content !== false) {
            $json = json_decode($content, true);
            if ($json !== null) {
                // Extraer algunas claves relevantes
                if (isset($json['resumen'])) {
                    $file_info['preview'] = $json['resumen'];
                } elseif (isset($json['organismo'])) {
                    $file_info['preview'] = ['organismo' => $json['organismo']];
                }
            }
        }
    }

    $results[] = $file_info;
}

// =============================================================================
// ORDENAR POR FECHA DE MODIFICACIÓN (MÁS RECIENTE PRIMERO)
// =============================================================================

usort($results, function($a, $b) {
    return $b['timestamp'] <=> $a['timestamp'];
});

// Remover timestamp del resultado final (solo para ordenar)
foreach ($results as &$result) {
    unset($result['timestamp']);
}

// =============================================================================
// RETORNAR RESULTADOS
// =============================================================================

debug_log("Found " . count($results) . " $type files");

json_success([
    'type' => $type,
    'count' => count($results),
    'results' => $results
]);

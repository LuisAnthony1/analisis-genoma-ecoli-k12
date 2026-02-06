<?php
/**
 * Endpoint: Buscar genomas en NCBI
 *
 * GET /backend/api/search_ncbi.php?query=escherichia&limit=20
 *
 * Parámetros:
 * - query (string): Término de búsqueda (organismo, accession, etc.)
 * - limit (int, opcional): Número máximo de resultados (default: 20, max: 50)
 *
 * Retorna:
 * - success: true/false
 * - query: término buscado
 * - count: número de resultados
 * - results: array de genomas encontrados
 */

require_once 'config.php';

// =============================================================================
// VALIDACIÓN DE PARÁMETROS
// =============================================================================

$query = isset($_GET['query']) ? sanitize_input($_GET['query']) : '';
$limit = isset($_GET['limit']) ? intval($_GET['limit']) : 20;

// Validar query
if (empty($query)) {
    json_error('El parámetro "query" es requerido', 400);
}

if (strlen($query) < 3) {
    json_error('La búsqueda debe tener al menos 3 caracteres', 400);
}

if (strlen($query) > 200) {
    json_error('La búsqueda no puede exceder 200 caracteres', 400);
}

// Validar limit
if ($limit < 1) {
    $limit = 1;
}
if ($limit > 50) {
    $limit = 50;
}

debug_log("Searching NCBI: query='$query', limit=$limit");

// =============================================================================
// EJECUTAR BÚSQUEDA
// =============================================================================

$result = execute_python('buscar_ncbi.py', [$query, $limit], 15);

// Verificar si hubo error en la ejecución
if (!$result['success']) {
    json_error(
        'Error al buscar en NCBI. Verifica tu conexión a internet.',
        500,
        $result['output']
    );
}

// Parsear el JSON retornado por Python
$data = json_decode($result['output'], true);

if ($data === null) {
    json_error(
        'Error al procesar respuesta del script Python',
        500,
        $result['output']
    );
}

// Verificar si el script Python retornó error
if (!$data['success']) {
    json_error(
        $data['error'] ?? 'Error desconocido al buscar en NCBI',
        500,
        $data
    );
}

// =============================================================================
// RETORNAR RESULTADOS
// =============================================================================

debug_log("Search completed: {$data['count']} results found");

json_success($data);

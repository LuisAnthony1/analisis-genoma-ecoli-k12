<?php
/**
 * Endpoint: Ejecutar script de análisis Python
 *
 * POST /backend/api/run_analysis.php
 *
 * Body JSON:
 * {
 *   "script": "analisis_genes",
 *   "organism": "ecoli_k12" (opcional, para scripts que lo requieren)
 * }
 *
 * Scripts permitidos:
 * - analisis_genes (requiere organism: 1=ecoli, 2=salmonella)
 * - analisis_codones (solo E. coli)
 * - comparar_genomas (requiere ambos genomas descargados)
 * - visualizaciones (requiere análisis previos)
 * - analisis_distancias_intergenicas (solo E. coli)
 *
 * Retorna:
 * - success: true/false
 * - output: salida del script
 * - execution_time: tiempo de ejecución en segundos
 */

require_once 'config.php';

// =============================================================================
// VALIDACIÓN DE PARÁMETROS
// =============================================================================

// Leer body JSON
$input_raw = file_get_contents('php://input');
$input = json_decode($input_raw, true);

if ($input === null) {
    json_error('Body JSON inválido', 400);
}

$script = isset($input['script']) ? sanitize_input($input['script']) : '';
$organism = isset($input['organism']) ? sanitize_input($input['organism']) : null;

// Scripts permitidos (whitelist de seguridad)
$allowed_scripts = [
    'analisis_genes' => [
        'file' => 'analisis_genes.py',
        'requires_organism' => true,
        'timeout' => 60
    ],
    'analisis_codones' => [
        'file' => 'analisis_codones.py',
        'requires_organism' => false,
        'timeout' => 30
    ],
    'comparar_genomas' => [
        'file' => 'comparar_genomas.py',
        'requires_organism' => false,
        'timeout' => 90
    ],
    'visualizaciones' => [
        'file' => 'visualizaciones.py',
        'requires_organism' => false,
        'timeout' => 120
    ],
    'analisis_distancias_intergenicas' => [
        'file' => 'analisis_distancias_intergenicas.py',
        'requires_organism' => false,
        'timeout' => 30
    ]
];

// Validar script
if (empty($script)) {
    json_error('El parámetro "script" es requerido', 400);
}

if (!isset($allowed_scripts[$script])) {
    json_error(
        'Script no permitido. Scripts válidos: ' . implode(', ', array_keys($allowed_scripts)),
        400
    );
}

$script_config = $allowed_scripts[$script];

// Mapear organism a número (según estructura en analisis_genes.py)
$organism_map = [
    'ecoli_k12' => '1',
    'salmonella_lt2' => '2'
];

$args = [];

// Validar organism si es requerido
if ($script_config['requires_organism']) {
    if (empty($organism)) {
        json_error('Este script requiere el parámetro "organism" (ecoli_k12 o salmonella_lt2)', 400);
    }

    if (!isset($organism_map[$organism])) {
        json_error('Organismo inválido. Valores válidos: ecoli_k12, salmonella_lt2', 400);
    }

    $args[] = $organism_map[$organism];
}

debug_log("Running analysis: $script with organism: " . ($organism ?? 'none'));

// =============================================================================
// EJECUTAR SCRIPT
// =============================================================================

$start_time = microtime(true);

$result = execute_python($script_config['file'], $args, $script_config['timeout']);

$execution_time = round(microtime(true) - $start_time, 2);

// =============================================================================
// RETORNAR RESULTADOS
// =============================================================================

// Nota: Algunos scripts pueden tener return_code != 0 pero completarse parcialmente
// No usamos $result['success'] como único criterio

debug_log("Analysis completed in {$execution_time}s - Return code: {$result['return_code']}");

json_success([
    'output' => $result['output'],
    'return_code' => $result['return_code'],
    'execution_time' => $execution_time,
    'script' => $script
]);

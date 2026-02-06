<?php
/**
 * Endpoint: Descargar genoma desde NCBI
 *
 * POST /backend/api/download_genome.php
 *
 * Body JSON:
 * {
 *   "accession_id": "NC_000913.3",
 *   "organism_name": "Escherichia coli K-12"
 * }
 *
 * Retorna:
 * - success: true/false
 * - message: mensaje descriptivo
 * - files: rutas de archivos generados (.gb, .fasta, metadata.json)
 * - metadata: información del genoma descargado
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

$accession_id = isset($input['accession_id']) ? sanitize_input($input['accession_id']) : '';
$organism_name = isset($input['organism_name']) ? sanitize_organism_name($input['organism_name']) : '';

// Validar accession ID
if (empty($accession_id)) {
    json_error('El parámetro "accession_id" es requerido', 400);
}

if (!validate_accession_id($accession_id)) {
    json_error('Formato de accession ID inválido. Ejemplo válido: NC_000913.3', 400);
}

// Validar organism name
if (empty($organism_name)) {
    json_error('El parámetro "organism_name" es requerido', 400);
}

if (strlen($organism_name) < 3 || strlen($organism_name) > 200) {
    json_error('El nombre del organismo debe tener entre 3 y 200 caracteres', 400);
}

debug_log("Downloading genome: $accession_id - $organism_name");

// =============================================================================
// EJECUTAR DESCARGA
// =============================================================================

// Timeout de 5 minutos (300s) para descargas de NCBI
$result = execute_python('descargar_genoma_api.py', [$accession_id, $organism_name], 300);

// Verificar si hubo error en la ejecución
if (!$result['success']) {
    json_error(
        'Error al descargar genoma desde NCBI',
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
        $data['error'] ?? 'Error desconocido al descargar genoma',
        500,
        $data
    );
}

// =============================================================================
// RETORNAR RESULTADOS
// =============================================================================

debug_log("Download completed: {$data['files']['genbank']}");

json_success($data, "Genoma descargado exitosamente: $organism_name");

<?php
/**
 * Endpoint: Listar genomas descargados
 *
 * GET /backend/api/list_genomes.php
 *
 * Retorna:
 * - success: true/false
 * - count: número de genomas descargados
 * - genomes: array de genomas con metadata
 */

require_once 'config.php';

// =============================================================================
// ESCANEAR ARCHIVOS GENBANK
// =============================================================================

$genomes = [];

// Buscar todos los archivos .gb en crudo/
$gb_files = glob(CRUDO_DIR . DIRECTORY_SEPARATOR . '*.gb');

if ($gb_files === false) {
    json_error('Error al escanear directorio de genomas', 500);
}

debug_log("Found " . count($gb_files) . " GenBank files");

// =============================================================================
// EXTRAER METADATA DE CADA GENOMA
// =============================================================================

foreach ($gb_files as $gb_file) {
    $basename = basename($gb_file, '.gb');
    $fasta_file = CRUDO_DIR . DIRECTORY_SEPARATOR . $basename . '.fasta';

    // Buscar archivo de metadata
    // Formato: metadata_ecoli.json o metadata_ecoli_k12.json
    $metadata_file = CRUDO_DIR . DIRECTORY_SEPARATOR . 'metadata_' . $basename . '.json';

    // Información básica del archivo
    $genome_info = [
        'filename' => basename($gb_file),
        'basename' => $basename,
        'size_mb' => round(filesize($gb_file) / 1048576, 2),
        'modified' => date('Y-m-d H:i:s', filemtime($gb_file)),
        'has_fasta' => file_exists($fasta_file),
        'has_metadata' => file_exists($metadata_file)
    ];

    // Leer metadata si existe
    if (file_exists($metadata_file)) {
        $metadata = json_decode(file_get_contents($metadata_file), true);

        if ($metadata !== null) {
            $genome_info['metadata'] = $metadata;

            // Extraer campos relevantes para mostrar en UI
            if (isset($metadata['genoma'])) {
                $genome_info['organism'] = $metadata['genoma']['organismo'] ?? 'Desconocido';
                $genome_info['accession_id'] = $metadata['id_acceso'] ?? 'N/A';
                $genome_info['length'] = $metadata['genoma']['longitud'] ?? 0;
                $genome_info['length_mb'] = round(($metadata['genoma']['longitud'] ?? 0) / 1e6, 2);
                $genome_info['download_date'] = $metadata['fecha_descarga'] ?? 'N/A';
                $genome_info['num_features'] = $metadata['genoma']['num_features'] ?? 0;
                $genome_info['description'] = $metadata['genoma']['descripcion'] ?? '';
            } else {
                // Metadata antigua sin estructura 'genoma'
                $genome_info['organism'] = $metadata['organismo_config'] ?? 'Desconocido';
                $genome_info['accession_id'] = $metadata['id_acceso'] ?? 'N/A';
                $genome_info['download_date'] = $metadata['fecha_descarga'] ?? 'N/A';
            }
        }
    } else {
        // No hay metadata, extraer info básica del nombre de archivo
        $genome_info['organism'] = ucwords(str_replace(['_', '-'], ' ', $basename));
        $genome_info['accession_id'] = 'N/A';
        $genome_info['download_date'] = date('Y-m-d H:i:s', filemtime($gb_file));
    }

    $genomes[] = $genome_info;
}

// =============================================================================
// ORDENAR POR FECHA DE DESCARGA (MÁS RECIENTE PRIMERO)
// =============================================================================

usort($genomes, function($a, $b) {
    return strtotime($b['download_date']) <=> strtotime($a['download_date']);
});

// =============================================================================
// RETORNAR RESULTADOS
// =============================================================================

debug_log("Returning " . count($genomes) . " genomes");

json_success([
    'count' => count($genomes),
    'genomes' => $genomes
]);

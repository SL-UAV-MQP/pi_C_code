/**
 * @file cell_tower_db.h
 * @brief Cell tower database for CMRCM Field area
 *
 * Maps Physical Cell IDs to known tower locations for direction finding.
 */

#ifndef CELL_TOWER_DB_H
#define CELL_TOWER_DB_H

#include <stdbool.h>
#include <stdint.h>

/* ============================================================================
 * Configuration Constants
 * ============================================================================ */

/** Maximum number of cell IDs per tower */
#define MAX_CELL_IDS_PER_TOWER 16

/** Maximum number of frequency bands per tower */
#define MAX_FREQ_BANDS 8

/** Maximum tower name length */
#define MAX_TOWER_NAME_LENGTH 64

/** Maximum carrier name length */
#define MAX_CARRIER_NAME_LENGTH 32

/** Maximum number of towers in database */
#define MAX_TOWERS 100

/* ============================================================================
 * Data Structures
 * ============================================================================ */

/**
 * @brief Cell tower information
 */
typedef struct {
    char name[MAX_TOWER_NAME_LENGTH];              /**< Tower name/identifier */
    double latitude;                                /**< Latitude in degrees */
    double longitude;                               /**< Longitude in degrees */
    int cell_ids[MAX_CELL_IDS_PER_TOWER];          /**< Physical Cell IDs */
    int num_cell_ids;                              /**< Number of Cell IDs */
    char carrier[MAX_CARRIER_NAME_LENGTH];         /**< Carrier name */
    int frequency_bands[MAX_FREQ_BANDS];           /**< Frequency bands (MHz) */
    int num_freq_bands;                            /**< Number of freq bands */
    double distance_from_cmrcm;                    /**< Distance in meters */
    double azimuth_from_cmrcm;                     /**< Azimuth in degrees */
} cell_tower_t;

/**
 * @brief Cell tower database
 */
typedef struct {
    cell_tower_t towers[MAX_TOWERS];               /**< Tower array */
    int num_towers;                                /**< Number of towers */

    // CMRCM Field reference coordinates
    double cmrcm_lat;                              /**< Reference latitude */
    double cmrcm_lon;                              /**< Reference longitude */
} cell_tower_database_t;

/* ============================================================================
 * Function Declarations
 * ============================================================================ */

/**
 * @brief Initialize cell tower database
 *
 * Sets CMRCM Field coordinates and initializes known towers.
 *
 * @param db Pointer to database structure
 * @return 0 on success, -1 on failure
 */
int cell_tower_db_init(cell_tower_database_t* db);

/**
 * @brief Calculate azimuth from CMRCM to a location
 *
 * Calculates bearing in degrees (0-360) from North.
 *
 * @param db Pointer to database structure
 * @param lat Target latitude in degrees
 * @param lon Target longitude in degrees
 * @return Azimuth in degrees (0-360)
 */
double cell_tower_db_calculate_azimuth(
    const cell_tower_database_t* db,
    double lat,
    double lon
);

/**
 * @brief Add tower to database
 *
 * @param db Pointer to database structure
 * @param tower Pointer to tower to add (copied into database)
 * @return 0 on success, -1 if database full
 */
int cell_tower_db_add_tower(
    cell_tower_database_t* db,
    const cell_tower_t* tower
);

/**
 * @brief Get tower by name
 *
 * @param db Pointer to database structure
 * @param name Tower name to search for
 * @return Pointer to tower or NULL if not found
 */
const cell_tower_t* cell_tower_db_get_by_name(
    const cell_tower_database_t* db,
    const char* name
);

/**
 * @brief Get tower by Physical Cell ID
 *
 * @param db Pointer to database structure
 * @param cell_id Physical Cell ID to search for
 * @return Pointer to tower or NULL if not found
 */
const cell_tower_t* cell_tower_db_get_by_cell_id(
    const cell_tower_database_t* db,
    int cell_id
);

/**
 * @brief Get tower by azimuth angle
 *
 * Finds tower within tolerance of specified azimuth.
 *
 * @param db Pointer to database structure
 * @param azimuth Azimuth in degrees (0-360)
 * @param tolerance Tolerance in degrees (default: 10.0)
 * @return Pointer to tower or NULL if not found
 */
const cell_tower_t* cell_tower_db_get_by_azimuth(
    const cell_tower_database_t* db,
    double azimuth,
    double tolerance
);

/**
 * @brief Get all towers in database
 *
 * @param db Pointer to database structure
 * @param towers Output array for tower pointers
 * @param max_towers Maximum towers to return
 * @return Number of towers returned
 */
int cell_tower_db_get_all_towers(
    const cell_tower_database_t* db,
    const cell_tower_t** towers,
    int max_towers
);

/**
 * @brief Get expected azimuth angles for all towers
 *
 * @param db Pointer to database structure
 * @param azimuths Output array for azimuths
 * @param max_azimuths Maximum azimuths to return
 * @return Number of azimuths returned
 */
int cell_tower_db_get_expected_azimuths(
    const cell_tower_database_t* db,
    double* azimuths,
    int max_azimuths
);

/**
 * @brief Associate a Cell ID with a tower
 *
 * Learns Cell ID mapping from field measurements.
 *
 * @param db Pointer to database structure
 * @param tower_name Tower name
 * @param cell_id Physical Cell ID to associate
 * @return 0 on success, -1 if tower not found or array full
 */
int cell_tower_db_associate_cell_id(
    cell_tower_database_t* db,
    const char* tower_name,
    int cell_id
);

/**
 * @brief Print database summary
 *
 * @param db Pointer to database structure
 */
void cell_tower_db_print_summary(const cell_tower_database_t* db);

#endif /* CELL_TOWER_DB_H */

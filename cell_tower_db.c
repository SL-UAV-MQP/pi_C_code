/**
 * @file cell_tower_db.c
 * @brief Cell tower database for CMRCM Field area
 *
 * Maps Physical Cell IDs to known tower locations for direction finding.
 */

#include "cell_tower_db.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* CMRCM Field reference coordinates */
#define CMRCM_LAT 42.307102
#define CMRCM_LON -71.612663

/* ============================================================================
 * Internal Helper Functions
 * ============================================================================ */

/**
 * @brief Calculate azimuth from CMRCM to a location
 * calculate_azimuth
 */
static double calculate_azimuth_internal(double lat1, double lon1, double lat2, double lon2) {
    // Convert to radians
    double lat1_rad = lat1 * M_PI / 180.0;
    double lon1_rad = lon1 * M_PI / 180.0;
    double lat2_rad = lat2 * M_PI / 180.0;
    double lon2_rad = lon2 * M_PI / 180.0;

    double dlon = lon2_rad - lon1_rad;

    double x = sin(dlon) * cos(lat2_rad);
    double y = cos(lat1_rad) * sin(lat2_rad) - sin(lat1_rad) * cos(lat2_rad) * cos(dlon);

    double azimuth = atan2(x, y) * 180.0 / M_PI;

    // Normalize to 0-360
    if (azimuth < 0) {
        azimuth += 360.0;
    }

    return azimuth;
}

/**
 * @brief Initialize known towers near CMRCM Field
 */
static void init_cmrcm_towers(cell_tower_database_t* db) {
    cell_tower_t tower;

    // Tower 1: WQNR305/WQNX339 - Closest tower (650m north)
    strcpy(tower.name, "WQNR305_WQNX339");
    tower.latitude = 42.312667;
    tower.longitude = -71.616139;
    tower.num_cell_ids = 0;
    strcpy(tower.carrier, "Unknown");
    tower.frequency_bands[0] = 700;
    tower.frequency_bands[1] = 850;
    tower.frequency_bands[2] = 1900;
    tower.frequency_bands[3] = 2100;
    tower.num_freq_bands = 4;
    tower.distance_from_cmrcm = 650.0;
    tower.azimuth_from_cmrcm = calculate_azimuth_internal(CMRCM_LAT, CMRCM_LON,
                                                          tower.latitude, tower.longitude);
    cell_tower_db_add_tower(db, &tower);

    // Tower 2: WQQL972/WQRV287 - Juniper Hill Golf Course (1.7km NW)
    strcpy(tower.name, "WQQL972_WQRV287");
    tower.latitude = 42.311194;
    tower.longitude = -71.630917;
    tower.num_cell_ids = 0;
    strcpy(tower.carrier, "Unknown");
    tower.frequency_bands[0] = 700;
    tower.frequency_bands[1] = 850;
    tower.frequency_bands[2] = 1900;
    tower.num_freq_bands = 3;
    tower.distance_from_cmrcm = 1700.0;
    tower.azimuth_from_cmrcm = calculate_azimuth_internal(CMRCM_LAT, CMRCM_LON,
                                                          tower.latitude, tower.longitude);
    cell_tower_db_add_tower(db, &tower);

    // Tower 3: Chauncy Hill (1.9km SE)
    strcpy(tower.name, "Chauncy_Hill");
    tower.latitude = 42.290389;
    tower.longitude = -71.615972;
    tower.num_cell_ids = 0;
    strcpy(tower.carrier, "Unknown");
    tower.frequency_bands[0] = 700;
    tower.frequency_bands[1] = 850;
    tower.frequency_bands[2] = 1900;
    tower.frequency_bands[3] = 2100;
    tower.num_freq_bands = 4;
    tower.distance_from_cmrcm = 1900.0;
    tower.azimuth_from_cmrcm = calculate_azimuth_internal(CMRCM_LAT, CMRCM_LON,
                                                          tower.latitude, tower.longitude);
    cell_tower_db_add_tower(db, &tower);
}

/* ============================================================================
 * Public API Implementation
 * ============================================================================ */

int cell_tower_db_init(cell_tower_database_t* db) {
    if (!db) {
        return -1;
    }

    db->num_towers = 0;
    db->cmrcm_lat = CMRCM_LAT;
    db->cmrcm_lon = CMRCM_LON;

    // Initialize known towers
    init_cmrcm_towers(db);

    return 0;
}

double cell_tower_db_calculate_azimuth(
    const cell_tower_database_t* db,
    double lat,
    double lon
) {
    // calculate_azimuth
    if (!db) {
        return 0.0;
    }

    return calculate_azimuth_internal(db->cmrcm_lat, db->cmrcm_lon, lat, lon);
}

int cell_tower_db_add_tower(
    cell_tower_database_t* db,
    const cell_tower_t* tower
) {
    // add_tower
    if (!db || !tower || db->num_towers >= MAX_TOWERS) {
        return -1;
    }

    db->towers[db->num_towers] = *tower;
    db->num_towers++;

    return 0;
}

const cell_tower_t* cell_tower_db_get_by_name(
    const cell_tower_database_t* db,
    const char* name
) {
    // get_tower_by_name
    if (!db || !name) {
        return NULL;
    }

    for (int i = 0; i < db->num_towers; i++) {
        if (strcmp(db->towers[i].name, name) == 0) {
            return &db->towers[i];
        }
    }

    return NULL;
}

const cell_tower_t* cell_tower_db_get_by_cell_id(
    const cell_tower_database_t* db,
    int cell_id
) {
    // get_tower_by_cell_id
    if (!db) {
        return NULL;
    }

    for (int i = 0; i < db->num_towers; i++) {
        for (int j = 0; j < db->towers[i].num_cell_ids; j++) {
            if (db->towers[i].cell_ids[j] == cell_id) {
                return &db->towers[i];
            }
        }
    }

    return NULL;
}

const cell_tower_t* cell_tower_db_get_by_azimuth(
    const cell_tower_database_t* db,
    double azimuth,
    double tolerance
) {
    // get_tower_by_azimuth
    if (!db) {
        return NULL;
    }

    for (int i = 0; i < db->num_towers; i++) {
        double angle_diff = fabs(db->towers[i].azimuth_from_cmrcm - azimuth);

        // Handle 360/0 degree wrap-around
        if (angle_diff > 180.0) {
            angle_diff = 360.0 - angle_diff;
        }

        if (angle_diff <= tolerance) {
            return &db->towers[i];
        }
    }

    return NULL;
}

int cell_tower_db_get_all_towers(
    const cell_tower_database_t* db,
    const cell_tower_t** towers,
    int max_towers
) {
    // get_all_towers
    if (!db || !towers) {
        return 0;
    }

    int count = (db->num_towers < max_towers) ? db->num_towers : max_towers;

    for (int i = 0; i < count; i++) {
        towers[i] = &db->towers[i];
    }

    return count;
}

int cell_tower_db_get_expected_azimuths(
    const cell_tower_database_t* db,
    double* azimuths,
    int max_azimuths
) {
    // get_expected_azimuths
    if (!db || !azimuths) {
        return 0;
    }

    int count = (db->num_towers < max_azimuths) ? db->num_towers : max_azimuths;

    for (int i = 0; i < count; i++) {
        azimuths[i] = db->towers[i].azimuth_from_cmrcm;
    }

    return count;
}

int cell_tower_db_associate_cell_id(
    cell_tower_database_t* db,
    const char* tower_name,
    int cell_id
) {
    // associate_cell_id
    if (!db || !tower_name) {
        return -1;
    }

    for (int i = 0; i < db->num_towers; i++) {
        if (strcmp(db->towers[i].name, tower_name) == 0) {
            // Check if cell_id already exists
            for (int j = 0; j < db->towers[i].num_cell_ids; j++) {
                if (db->towers[i].cell_ids[j] == cell_id) {
                    return 0;  // Already associated
                }
            }

            // Add new cell_id if space available
            if (db->towers[i].num_cell_ids < MAX_CELL_IDS_PER_TOWER) {
                db->towers[i].cell_ids[db->towers[i].num_cell_ids] = cell_id;
                db->towers[i].num_cell_ids++;
                return 0;
            }

            return -1;  // Array full
        }
    }

    return -1;  // Tower not found
}

void cell_tower_db_print_summary(const cell_tower_database_t* db) {
    // print_summary
    if (!db) {
        return;
    }

    printf("\n=== CMRCM Field Cell Tower Database ===\n");
    printf("Reference: %.6f, %.6f\n\n", db->cmrcm_lat, db->cmrcm_lon);

    // Sort by distance (simple bubble sort for small arrays)
    int* indices = (int*)malloc(db->num_towers * sizeof(int));
    for (int i = 0; i < db->num_towers; i++) {
        indices[i] = i;
    }

    for (int i = 0; i < db->num_towers - 1; i++) {
        for (int j = 0; j < db->num_towers - i - 1; j++) {
            if (db->towers[indices[j]].distance_from_cmrcm >
                db->towers[indices[j+1]].distance_from_cmrcm) {
                int temp = indices[j];
                indices[j] = indices[j+1];
                indices[j+1] = temp;
            }
        }
    }

    for (int i = 0; i < db->num_towers; i++) {
        const cell_tower_t* tower = &db->towers[indices[i]];

        printf("%s:\n", tower->name);
        printf("  Location: %.6f, %.6f\n", tower->latitude, tower->longitude);
        printf("  Distance: %.0fm\n", tower->distance_from_cmrcm);
        printf("  Azimuth: %.1f°\n", tower->azimuth_from_cmrcm);
        printf("  Carrier: %s\n", tower->carrier);

        printf("  Cell IDs: ");
        if (tower->num_cell_ids > 0) {
            printf("[");
            for (int j = 0; j < tower->num_cell_ids; j++) {
                printf("%d", tower->cell_ids[j]);
                if (j < tower->num_cell_ids - 1) {
                    printf(", ");
                }
            }
            printf("]\n");
        } else {
            printf("Unknown\n");
        }

        printf("\n");
    }

    free(indices);
}

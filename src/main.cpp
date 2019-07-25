#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include "getopt.h"

int main( int argc, char *argv[] ) {
    int opt;
    uint particles, iterations;
    char* matrix_file = 0;
    char* mutations_file = 0;
    std::vector<std::vector<int>> matrix;

    uint cells, mutations;
    std::vector<std::string> mutation_names;

    while ( ( opt = getopt( argc, argv, "p:i:f:m:" ) ) != -1 ) {
        switch ( opt ) {
        case 'p': // particles
            particles = atoi( optarg );
            break;
        case 'i': // iterations
            iterations = atoi( optarg );
            break;
        case 'f': // matrix input file
            matrix_file = optarg;
            break;
        case 'm': // mutations input file
            mutations_file = optarg;
            break;
        case ':':
            switch( optopt ) {
                case 'm':
                    mutations_file = 0;
                    break;
                default:
                    fprintf( stderr, "Option -%c is missing a required argument\n", optopt );
                    exit( EXIT_FAILURE );
                    break;
            }
            break;
        default:
            fprintf( stderr, "Usage: %s -p particles -i iterations -i matrix_file -m mutation_names\n", argv[0] );
            exit( EXIT_FAILURE );
            break;
        }
    }

    // reading from the raw data file
    std::ifstream infile( matrix_file );
    if ( !infile ) {
        printf( "Filename %s not found!", matrix_file );
        exit( EXIT_FAILURE );
    }

    char value;
    std::string line;

    while ( std::getline( infile, line ) ) {

        std::istringstream linestream(line);
        std::vector<int> value_line;
        while ( linestream >> value ) {
            value_line.push_back( value );
        }
        matrix.push_back( value_line );
    }

    cells = matrix.size();

    if ( cells < 1 ) {
        printf( "Error: no cells were found in the file." );
        exit( EXIT_FAILURE );
    }

    mutations = matrix.at(0).size();

    // Mutations
    if ( mutations_file == 0 ) {

        for ( uint i = 1; i <= mutations; i++ ) {
            mutation_names.push_back( std::to_string( i ) );
        }

    } else {
        // reading from the raw data file
        std::ifstream infile( mutations_file );
        if ( !infile ) {
            printf( "Filename %s not found!", mutations_file );
            exit( EXIT_FAILURE );
        }

        std::string line;

        while ( std::getline( infile, line ) ) {
            mutation_names.push_back( line );
        }

        if ( mutations < mutation_names.size() ) {
            fprintf( stderr, "Mutations in the data file and mutation names file do not match!" );
            exit( EXIT_FAILURE );
        }
    }

    

    return 0;
    
}

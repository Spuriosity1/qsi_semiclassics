#ifndef PARALLEL_IO_HH
#define PARALLEL_IO_HH

#include <chrono>
#include <thread>
#include <memory>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <fstream>

template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    int size_s = std::snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    if( size_s <= 0 ){ throw std::runtime_error( "Error during formatting." ); }
    auto size = static_cast<size_t>( size_s );
    std::unique_ptr<char[]> buf( new char[ size ] );
    std::snprintf( buf.get(), size, format.c_str(), args ... );
    return std::string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}


/**
 * @brief Writes "what" to file "filename" without colliding
 * 
 * Will check if there is a file by the name [filename]+[suffix] (default suffix is .lock)
 * if the lockfile exists, it will wait a random amount of time between 0 and 10ms
 * (to prevent cycling) and try again.
 * 
 * Function will give up if success does not happen within nattempt tries.
 * 
 * 
 * @param filename 
 * @param what 
 * @param suffix 
 * @param nattempt 
 */
void print_without_collision(const std::filesystem::path& filename, const std::string& what,
const char* suffix=".lock", unsigned long nattempt=100){
    //check if lockfile exists
    std::filesystem::path lockfile = filename;
    lockfile += suffix;

    for (unsigned long ncall=0; ncall < nattempt; ncall++){
        if (std::filesystem::exists(lockfile)){
            // bad news
            if (ncall > nattempt){
                throw std::runtime_error("Exceeded maximum number of attempts to save data, \
                maybe a process died leaving a lockfile behind?");
            } else {
                printf("Waiting my turn (%4ld/%4ld)\n", ncall, nattempt);
                std::this_thread::sleep_for(std::chrono::microseconds(rand()%10000));
            }
        }
    }

    std::ofstream lock_of(lockfile);
    lock_of << "ncall\n";
    std::ofstream of(filename, std::ios_base::app);
    of << what;
    of.close();
    lock_of.close();
    // allow other instances to append by deleting the lockfile
    std::filesystem::remove(lockfile);
}

#endif
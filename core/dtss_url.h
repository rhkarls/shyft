#pragma once
#include <algorithm>
#include <unordered_map>
#include <string>

namespace shyft {
namespace dtss {

// TODO: inline when vs implements P0386R2: Inline variables
extern std::string shyft_prefix;//="shyft://";  ///< marks all internal handled urls


/**construct a shyft-url from container and ts-name */
inline std::string shyft_url(const std::string& container, const std::string& ts_name) {
    return shyft_prefix + container + "/" + ts_name;
}


/** match & extract fast the following 'shyft://<container>/'
 * \param url like pattern above
 * \return <container> or empty string if no match
 */
inline std::string extract_shyft_url_container(const std::string& url) {
    if ( (url.size() < shyft_prefix.size() + 2) || !std::equal(begin(shyft_prefix), end(shyft_prefix), begin(url)) )
        return std::string{};
    // path after container?
    auto ce = url.find_first_of('/', shyft_prefix.size());  // container end
    if ( ce == std::string::npos )
        return std::string{};
    return url.substr(shyft_prefix.size(), ce - shyft_prefix.size());
}

/** Extract any query parameters from a url.
 *
 * The implementation ignores data until the first `?` character, afterwhich it parses the rest of
 * the url as a query string.
 *
 * The query string is assumed to be on the format `?key1=value1&key2=value2&key3=&key4=value4`.
 * This will be parsed into a map with four keys: `key1` through `key4`, where `key3` have a
 * empty string value, while the rest have respectivly values `value1`, `value2`, and `value4`.
 *
 * If the url does not have any query parameters an empty mapping wil be returned.
 *
 * \param url String url to parse.
 * \return A map-type with the key-values from the url query string.
 */
inline std::unordered_map<std::string, std::string> extract_query_parameters(const std::string & url) {
    
    using map_t = std::unordered_map<std::string, std::string>;
    
    // locate query string if present
    auto it = std::find(url.cbegin(), url.cend(), '?');
    if ( it == url.cend() )
        return map_t{};
    std::advance(it, 1);  // skip to first character of the query string
    if ( it == url.cend() )
        return map_t{};

    // parse key=value pairs
    map_t queries{};
    std::string key_str{};
    std::string val_str{};
    while ( it != url.cend() ) {
        // locate key
        auto end = std::find(it, url.cend(), '=');
        if ( end == url.cend() ) {
            // no = after key, and at the end of the string -> do not add to map, then return
            break;
        }
        key_str.assign(it, end);  // found key!

        // locate value
        it = end + 1;  // skip pask =
        end = std::find(it, url.cend(), '&');
        if ( std::distance(it, end) == 0 ) {
            // no value after = -> use empty string
            val_str.assign("");
        } else {
            // value after =
            val_str.assign(it, end);
        }

        // store and advance
        queries[key_str] = val_str;
        it = end;
        if ( it != url.cend() ) {
            std::advance(it, 1);
        }
    }
    return queries;
}

}
}
#pragma once

#include <algorithm>
#include <iterator>
#include <map>
#include <sstream>
#include <string>


namespace shyft {
namespace dtss {

template < std::size_t N >
inline constexpr std::size_t char_str_length(const char(&a)[N]) {
    return N-1;  // don't count \0 char
}


// TODO: inline when vs implements P0386R2: Inline variables
constexpr char shyft_prefix[] = "shyft://";  ///< marks all internal handled urls


/** Construct a shyft-url from a container and a ts-name. */
inline std::string shyft_url(const std::string & container, const std::string & ts_name) {
    return std::string{ shyft_prefix } + container + "/" + ts_name;
}
/** Construct a shyft-url from a container, a ts-name, and a collection of query flags. */
inline std::string shyft_url(const std::string & container, const std::string & ts_name, const std::map<std::string, std::string> & queries) {
    std::ostringstream str_s{ };
    str_s << "?";
    for ( auto it = queries.cbegin(); it != queries.cend(); ) {
        auto p = *it++;  // record current value, increment iterator
        str_s << p.first << '=' << p.second << (it != queries.cend() ? "&" : "");
    }
    return shyft_url(container, ts_name) + str_s.str();
}


/** match & extract fast the following 'shyft://<container>/'
 * \param url like pattern above
 * \return <container> or empty string if no match
 */
inline std::string extract_shyft_url_container(const std::string & url) {
    if ( (url.size() < char_str_length(shyft_prefix) + 2) || !std::equal(std::begin(shyft_prefix), std::prev(std::end(shyft_prefix)), begin(url)) )
        return std::string{};
    // path after container?
    auto ce = url.find_first_of('/', char_str_length(shyft_prefix));  // container end
    if ( ce == std::string::npos )
        return std::string{};
    return url.substr(char_str_length(shyft_prefix), ce - char_str_length(shyft_prefix));
}


inline std::string extract_shyft_url_path(const std::string & url) {
    if ( url.size() < (char_str_length(shyft_prefix) + 2) || ! std::equal(std::begin(shyft_prefix), std::prev(std::end(shyft_prefix)), begin(url)) )
        return std::string{};
    auto ce = url.find_first_of('/', char_str_length(shyft_prefix));  // container end
    if ( ce == url.npos )
        return std::string{};
    auto qs = url.find('?', ce);
    if ( qs == url.npos ) {
        return url.substr(ce + 1);
    } else {
        return url.substr(ce + 1, qs - ce - 1);
    }
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
inline std::map<std::string, std::string> extract_shyft_url_query_parameters(const std::string & url) {
    
    using map_t = std::map<std::string, std::string>;
    
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
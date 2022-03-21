// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************

#ifndef _AFLOWLIB_ENTRY_LOADER_H_
#define _AFLOWLIB_ENTRY_LOADER_H_


namespace aflowlib {
  class EntryLoader : public xStream {
    public:
      //NECESSARY PUBLIC CLASS METHODS - START
      //constructors - START
      EntryLoader(std::ostream& oss=cout);
      EntryLoader(std::ofstream& FileMESSAGE,std::ostream& oss=cout);
      EntryLoader(const EntryLoader& b);
      //constructors - STOP
      ~EntryLoader();
      const EntryLoader& operator=(const EntryLoader& other);
      void clear();
      //NECESSARY PUBLIC CLASS METHODS - STOP

      /// @brief available sources to load lib entries
      enum class Source {
        SQLITE,         ///< internal AFLUX database
        AFLUX,          ///< public AFLUX API
        FILESYSTEM,     ///< search by public alloy DB, loaded from internal filesystem
        FILESYSTEM_RAW, ///< internal filesystem
        RESTAPI,        ///< search by public alloy DB, loaded from public REST API
        RESTAPI_RAW,    ///< public REST API
        NONE,           ///< no source set
        FAILED          ///< setting a source failed
      };

      // Settings
      bool m_out_silent = false; ///< silents all output
      bool m_out_debug = false;  ///< verbose output (overwrites #m_out_silent)

      bool m_xstructure_relaxed= false;  ///< add the relaxed structure into the lib entry
      bool m_xstructure_original= false; ///< add the original structure into the lib entry

      /// alternative sources for the original structure if aflow.in failed (in order of priority)
      std::vector<std::string> m_xstructure_original_file_name = {"POSCAR.orig", "POSCAR.relax1"};
      /// sources for the relaxed structure (in order of priority)
      std::vector<std::string> m_xstructure_final_file_name = {"CONTCAR.relax2", "CONTCAR.relax", "POSCAR.static", "POSCAR.bands", "CONTCAR.static", "CONTCAR.bands"};
      
      std::string m_sqlite_file = DEFAULT_AFLOW_DB_FILE;                    ///< location of the internal AFLUX SQLITE DB file
      std::string m_sqlite_alloy_file = DEFAULT_ENTRY_LOADER_ALLOY_DB_FILE; ///< location of the public alloy SQLITE DB file
      std::string m_sqlite_collection = "WEB";                              ///< collection to use for queries
      
      std::string m_aflux_server = DEFAULT_ENTRY_LOADER_AFLUX_SERVER;       ///< server that provides the AFLUX API
      std::string m_aflux_path = DEFAULT_ENTRY_LOADER_AFLUX_PATH;           ///< base path including AFLUX API version
      std::string m_aflux_collection = "RAW";                               ///< collection to use for queries
      /// default AFLUX API directives to set format and disable paging
      std::map<std::string, std::string> m_aflux_directives {{"format", "aflow"}, {"paging", "0"}};

      std::string m_restapi_server = DEFAULT_ENTRY_LOADER_RESTAPI_SERVER;   ///< server that provides the AFLOW RESTAPI
      std::string m_restapi_path = DEFAULT_ENTRY_LOADER_RESTAPI_PATH;       ///< AFLOW RESTAPI base path
      std::string m_restapi_directives = DEFAULT_FILE_AFLOWLIB_ENTRY_OUT;   ///< directive to receive an entry (typically "?format=text" or "aflowlib.out")
      std::string m_restapi_listing_dirs = "?aflowlib_entries";             ///< directive to list sub entries
      std::string m_restapi_listing_files = "?files";                       ///< directive to list files available for an entry
      std::string m_restapi_collection = "WEB";                             ///< collection to use for queries

      std::string m_filesystem_outfile = DEFAULT_FILE_AFLOWLIB_ENTRY_OUT;   ///< file name of the AFLOW lib entry
      std::string m_filesystem_path = DEFAULT_ENTRY_LOADER_FS_PATH;         ///< base path to the internal filesystem
      std::string m_filesystem_collection = "RAW";                          ///< collection to use for queries

      // Attributes
      std::shared_ptr<aflowlib::AflowDB> m_sqlite_db_ptr;                   ///< pointer to an instance of the internal AFLUX SQLITE DB
      std::shared_ptr<aflowlib::AflowDB> m_sqlite_alloy_db_ptr;             ///< pointer to an instance of the public alloy SQLITE DB
      /// @brief set that contains all loaded AUIDs
      /// @note `std::set` is stored sorted and is therefore faster at finding entries compared to `std::vector`
      std::unordered_set<std::string> m_auid_list;

      // Data views
      /// @brief flat collection of loaded lib entries
      std::shared_ptr<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>> m_entries_flat;
      /// @brief tiered collection of loaded lib entries (number of elements -> alloy -> entries)
      /// @note creating a copy of this smart pointer is the most efficient way to use
      ///       the loaded entries after the EntryLoader class goes out of scope
      std::shared_ptr<std::map<short, std::map<std::string, std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>>> m_entries_layered_map;

      // REGEX expressions for quick finding/replacements in strings
      /// REGEX to find all chemical elements in a string
      const std::regex m_re_elements{"(A[cglmrstu]|B[aehikr]?|C[adeflmnorsu]?|D[bsy]|E[rsu]|F[elmr]?|G[ade]|H[efgos]?|I[nr]?|Kr?|L[airuv]|M[dgnot]|N[abdeiop]?|Os?|P[abdmortu]?|R[abefghnu]|S[bcegimnr]?|T[abcehilm]|U(u[opst])?|V|W|Xe|Yb?|Z[nr])"};
      /// REGEX to replace all pseudo potentials that contain uppercase letters from a string (could be miss taken for a chemical element)
      const std::regex m_re_ppclean{"(_GW|_AE|_200eV)"};
      /// @brief REGEX to help change a AURL into a file path
      /// @note the content of the group `((?:(?:LIB\d{1,})|(?:ICSD)))` can be used in the replacement  with `$1`;
      ///       the second group `(?:(?:RAW)|(?:LIB)|(?:WEB))` is there to select the full substring to be replaced
      const std::regex m_re_aurl2file{"((?:(?:LIB\\d{1,})|(?:ICSD)))_(?:(?:RAW)|(?:LIB)|(?:WEB))\\/"};

      // Generic entry loaders
      void loadAUID(std::string AUID);
      void loadAUID(const std::vector<std::string> &AUID);

      void loadAURL(std::string AURL);
      void loadAURL(const std::vector<std::string> &AURL);
      
      void loadAlloy(const std::string & alloy, bool recursive=true);
      void loadAlloy(const std::vector<std::string> &alloy, bool recursive=true);

      //TODO add loadICSD()

      // Direct entry loader
      void loadSqliteWhere(const std::string & where);
      void loadAFLUXQuery(const std::string & query);
      void loadAFLUXMatchbook(const std::map<std::string, std::string> & matchbook);
      void loadRestAPIQueries(const std::vector<std::string> & queries, bool full_url=false);
      void loadFiles(const std::vector<std::string> & files);
      void loadText(const std::vector<std::string> & raw_data_lines);

      // Source setter and getter
      bool setSource(EntryLoader::Source new_source);
      EntryLoader::Source getSource() const;

      // Getter for raw data
      std::vector<std::string> getRawSqliteWhere(const std::string & where) const;
      std::vector<std::string> getRawAFLUXMatchbook(const std::map<std::string, std::string> & matchbook);
      std::vector<std::string> getRawAFLUXQuery(const std::string & query);
      std::string getRawRestAPIQuery(const std::string &query, bool full_url=false);

      // xstructure loaders
      void addXstructure(aflowlib::_aflowlib_entry & entry, bool orig = false);
      bool loadXstructureFile(const aflowlib::_aflowlib_entry & entry, xstructure & new_structure, bool orig=false);
      bool loadXstructureLibEntry(const aflowlib::_aflowlib_entry & entry, xstructure & new_structure);
      bool loadXstructureAflowIn(const aflowlib::_aflowlib_entry & entry, xstructure & new_structure);

      // getter for views
      void getEntriesViewFlat(std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>> & result) const;
      void getEntriesViewTwoLayer(vector<vector<std::shared_ptr<aflowlib::_aflowlib_entry>>> & result) const;
      void getEntriesViewThreeLayer(std::vector<std::vector<std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>> & result) const;
      void getEntriesViewMap(std::map<short,std::map<std::string,std::vector<std::shared_ptr<aflowlib::_aflowlib_entry>>>> & result) const;

      // getter that copy data into new format
      void getEntriesFlat(std::vector<aflowlib::_aflowlib_entry> & result) const;
      void getEntriesTwoLayer(std::vector<std::vector<aflowlib::_aflowlib_entry>> & result) const;
      void getEntriesThreeLayer(std::vector<std::vector<vector<aflowlib::_aflowlib_entry>>> & result) const;

      
    private:
      //NECESSARY private CLASS METHODS - START
      void init();
      void copy(const EntryLoader& b);
      //NECESSARY END CLASS METHODS - END

      /// @brief Determines which source is used to load data
      /// @note should just be set by setSource()
      Source m_current_source = Source::NONE;
      bool m_filesystem_available = false; ///< helper to avoid repeatedly testing if the filesystem is available
      std::stringstream m_logger_message;  ///< reusable stringstream for logging
      bool m_out_super_silent = false;     ///< mute all output - used for private checks

      /// available ICSD symmetries to speed up the REST API searches
      const std::vector<std::string> m_icsd_symmetries = {"BCC/","BCT/","CUB/","FCC/","HEX/","MCL/","MCLC/",
                                                          "ORC/","ORCC/","ORCF/","ORCI/","RHL/","TET/","TRI/"};

      void selectSource();
      void listRestAPI(std::string url, std::vector<std::string> & result, bool directories=true);
      void getAlloyAUIDList(const std::vector<std::string> & alloy_list, std::vector<std::string> & auid_list);
      void loadAlloySearchFSR(const std::vector<std::string> & alloy_list, uint lib_max); // FILESYSTEM_RAW
      void loadAlloySearchRR(const std::vector<std::string> & alloy_list, uint lib_max); // RESTAPI_RAW
      bool cleanAUID(std::string & AUID);
      bool cleanAURL(std::string & AURL);
      std::string buildAFLUXQuery(const std::map<std::string, std::string> & matchbook) const;
      std::string extractAlloy(std::string name, char lib_type) const;
      std::string aflowin2poscar(std::string raw_in) const;

      // Logging helper
      void outInfo(const std::string & function_name);
      void outDebug(const std::string & function_name);
      void outError(const std::string & function_name, int line_number);
      void outHardError(const std::string & function_name, int line_number, int error_type);

  };
} // namespace aflowlib

#endif  // _AFLOWLIB_ENTRY_LOADER_H_

// ***************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2022           *
// *                                                                         *
// ***************************************************************************

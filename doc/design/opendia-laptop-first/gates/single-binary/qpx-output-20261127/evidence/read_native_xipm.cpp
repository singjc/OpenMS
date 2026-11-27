#include <OpenMS/FORMAT/XIPMParquetFile.h>
#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>

int main(int argc, char** argv)
{
  if (argc != 2) return 2;
  try
  {
    OpenMS::XIPMParquetFile reader(argv[1]);
    std::vector<OpenMS::XIPMParquetFile::XIPMPeakMap> maps;
    reader.load(maps);
    std::set<OpenMS::Int64> precursors;
    std::size_t points = 0, ms1 = 0, ms2 = 0;
    for (const auto& map : maps)
    {
      const auto n = map.mz.size();
      if (map.rt.size() != n || map.intensity.size() != n || map.ion_mobility.size() != n)
        throw std::runtime_error("Native XIPM reader returned unequal array lengths");
      if (!map.has_precursor_id || !map.has_precursor_decoy || map.precursor_decoy != 0)
        throw std::runtime_error("Missing precursor metadata or unexpected decoy");
      precursors.insert(map.precursor_id);
      if (map.ms_level == 1) ++ms1;
      else if (map.ms_level == 2) ++ms2;
      else throw std::runtime_error("Unexpected MS level");
      for (std::size_t i = 0; i < n; ++i)
      {
        if (!std::isfinite(map.mz[i]) || !std::isfinite(map.rt[i]) || !std::isfinite(map.intensity[i]) || map.intensity[i] < 0)
          throw std::runtime_error("Invalid measured point");
        if (i && map.rt[i] < map.rt[i - 1]) throw std::runtime_error("Nonmonotone map RT");
      }
      points += n;
    }
    if (maps.empty() || points == 0 || !ms1 || !ms2) throw std::runtime_error("Representative XIPM unexpectedly empty");
    std::cout << "NATIVE_XIPM_PASS maps=" << maps.size() << " precursors=" << precursors.size()
              << " ms1=" << ms1 << " ms2=" << ms2 << " points=" << points << '\n';
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}

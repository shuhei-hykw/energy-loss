// g4_table_generator
//
// Minimal Geant4 CLI that tabulates stopping power and CSDA range for one
// (particle, material) combination across a user-supplied energy grid.
// Uses G4EmCalculator backed by G4EmStandardPhysics_option4 (the same
// physics list NIST PSTAR/ASTAR are aligned with).
//
// Output: NIST-PSTAR-compatible CSV so the same Python
// RangeEnergyTable.from_nist_csv() parser can read it. The metadata
// header records Geant4 version, physics list, material density and
// composition, particle, and date.
//
// Usage:
//   g4_table_generator
//     --particle <name>       Geant4 particle name (e.g. "proton", "alpha")
//     --material <G4 NIST>    Geant4 NIST material name (e.g. "G4_PHOTO_EMULSION")
//     --emin <MeV>            Lower bound of the energy grid
//     --emax <MeV>            Upper bound of the energy grid
//     --n <int>               Number of points (>=2)
//     --grid <log|linear>     Default: log
//     --output <path>         CSV path
//
// All energies are MeV in/out; ranges are written in g/cm^2 to match the
// NIST PSTAR CSV format that the Python side already parses.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "G4Box.hh"
#include "G4EmCalculator.hh"
#include "G4EmParameters.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4Element.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysListFactory.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Version.hh"

namespace {

// Minimal modular physics list wrapping the EM Standard option4 physics
// constructor — same EM physics line PSTAR / ASTAR are aligned with.
class EmOption4PhysicsList : public G4VModularPhysicsList {
public:
  EmOption4PhysicsList() {
    RegisterPhysics(new G4EmStandardPhysics_option4());
  }
};

class Detector : public G4VUserDetectorConstruction {
public:
  explicit Detector(const std::string& mat) : mat_(mat) {}
  G4VPhysicalVolume* Construct() override {
    auto* nist = G4NistManager::Instance();
    G4Material* mat = nist->FindOrBuildMaterial(mat_);
    if (!mat) {
      std::cerr << "ERROR: unknown Geant4 NIST material: " << mat_ << "\n";
      std::exit(2);
    }
    // 1 m vacuum world hosting a single block of the target material.
    G4Material* vac = nist->FindOrBuildMaterial("G4_Galactic");
    auto* worldSolid = new G4Box("world", 1.0 * m, 1.0 * m, 1.0 * m);
    auto* worldLog = new G4LogicalVolume(worldSolid, vac, "world");
    auto* worldPhys = new G4PVPlacement(
      nullptr, {}, worldLog, "world", nullptr, false, 0);
    auto* sampleSolid = new G4Box("sample", 0.5 * m, 0.5 * m, 0.5 * m);
    auto* sampleLog = new G4LogicalVolume(sampleSolid, mat, "sample");
    new G4PVPlacement(
      nullptr, {}, sampleLog, "sample", worldLog, false, 0);
    material_ = mat;
    return worldPhys;
  }
  G4Material* material() const { return material_; }
private:
  std::string mat_;
  G4Material* material_ = nullptr;
};

struct Args {
  std::string particle = "proton";
  std::string material = "G4_PHOTO_EMULSION";
  double emin_mev = 1.0e-3;
  double emax_mev = 1.0e4;
  int n_points = 100;
  std::string grid = "log";
  std::string output = "g4_table.csv";
};

Args parse_args(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    std::string key = argv[i];
    auto need = [&](const char* k) {
      if (i + 1 >= argc) {
        std::cerr << "missing value for " << k << "\n";
        std::exit(2);
      }
      return std::string(argv[++i]);
    };
    if (key == "--particle") a.particle = need("--particle");
    else if (key == "--material") a.material = need("--material");
    else if (key == "--emin") a.emin_mev = std::stod(need("--emin"));
    else if (key == "--emax") a.emax_mev = std::stod(need("--emax"));
    else if (key == "--n") a.n_points = std::stoi(need("--n"));
    else if (key == "--grid") a.grid = need("--grid");
    else if (key == "--output") a.output = need("--output");
    else {
      std::cerr << "unknown argument: " << key << "\n";
      std::exit(2);
    }
  }
  if (a.n_points < 2) {
    std::cerr << "--n must be >= 2\n";
    std::exit(2);
  }
  if (!(a.emin_mev > 0.0 && a.emax_mev > a.emin_mev)) {
    std::cerr << "require 0 < emin < emax\n";
    std::exit(2);
  }
  if (a.grid != "log" && a.grid != "linear") {
    std::cerr << "--grid must be 'log' or 'linear'\n";
    std::exit(2);
  }
  return a;
}

std::vector<double> make_grid(const Args& a) {
  std::vector<double> ts(a.n_points);
  if (a.grid == "log") {
    const double lo = std::log(a.emin_mev);
    const double hi = std::log(a.emax_mev);
    for (int i = 0; i < a.n_points; ++i) {
      const double f = static_cast<double>(i) / (a.n_points - 1);
      ts[i] = std::exp(lo + f * (hi - lo));
    }
  } else {
    const double step = (a.emax_mev - a.emin_mev) / (a.n_points - 1);
    for (int i = 0; i < a.n_points; ++i) {
      ts[i] = a.emin_mev + i * step;
    }
  }
  return ts;
}

std::string clean_geant4_version() {
  // G4VERSION_NUMBER = major*100 + minor*10 + patch (e.g. 1141 -> 11.4.1).
  // We use this instead of the G4Version string macro because the
  // latter still carries SVN $Name:...$ keywords that look noisy in
  // the CSV provenance header.
  const int n = G4VERSION_NUMBER;
  const int major = n / 100;
  const int minor = (n / 10) % 10;
  const int patch = n % 10;
  std::ostringstream os;
  os << major << '.' << minor << '.' << patch;
  return os.str();
}

std::string today_iso() {
  std::time_t now = std::time(nullptr);
  std::tm tm{};
  localtime_r(&now, &tm);
  char buf[16];
  std::strftime(buf, sizeof(buf), "%Y-%m-%d", &tm);
  return buf;
}

std::string composition_string(const G4Material* mat) {
  std::ostringstream os;
  const auto* elems = mat->GetElementVector();
  const auto* fracs = mat->GetFractionVector();
  os << std::setprecision(6) << std::fixed;
  for (std::size_t i = 0; i < mat->GetNumberOfElements(); ++i) {
    if (i) os << ", ";
    os << (*elems)[i]->GetSymbol() << "=" << fracs[i];
  }
  return os.str();
}

} // namespace


int main(int argc, char** argv) {
  Args args = parse_args(argc, argv);

  // Run manager (serial — we don't need MT for table generation).
  auto* runManager = G4RunManagerFactory::CreateRunManager(
    G4RunManagerType::SerialOnly);

  // BuildCSDARange must be enabled before physics tables are built,
  // otherwise G4EmCalculator::GetCSDARange returns zero.
  G4EmParameters::Instance()->SetBuildCSDARange(true);

  auto* det = new Detector(args.material);
  runManager->SetUserInitialization(det);
  runManager->SetUserInitialization(new EmOption4PhysicsList());
  runManager->Initialize();
  // BeamOn(0) forces the physics tables to be built (cross sections,
  // dE/dx, range) for the registered particles in the registered
  // materials. Without it G4EmCalculator can crash on table lookups.
  runManager->BeamOn(0);

  G4Material* mat = det->material();
  if (!mat) {
    std::cerr << "ERROR: detector failed to construct material\n";
    return 3;
  }
  auto* part =
    G4ParticleTable::GetParticleTable()->FindParticle(args.particle);
  if (!part) {
    std::cerr << "ERROR: unknown Geant4 particle: " << args.particle << "\n";
    return 4;
  }

  G4EmCalculator calc;
  const auto ts_mev = make_grid(args);

  std::ofstream out(args.output);
  if (!out) {
    std::cerr << "ERROR: cannot open output: " << args.output << "\n";
    return 5;
  }

  const double rho_g_per_cm3 = mat->GetDensity() / (g / cm3);
  const double i_ev = mat->GetIonisation()->GetMeanExcitationEnergy() / eV;

  const std::string g4ver = clean_geant4_version();
  out << "# Geant4 stopping-power / CSDA-range table\n"
      << "# Particle: " << args.particle << "\n"
      << "# Source:   geant4://" << g4ver << "\n"
      << "# Geant4 version: " << g4ver << "\n"
      << "# Underlying physics: G4EmStandardPhysics_option4\n"
      << "# Density [g/cm^3]: " << std::fixed << std::setprecision(5)
      << rho_g_per_cm3 << "\n"
      << "# Mean excitation energy I [eV]: " << std::setprecision(3)
      << i_ev << "\n"
      << "# Composition (weight fraction): " << composition_string(mat)
      << "\n"
      << "# NIST table for " << mat->GetName() << "\n"
      << "# Fetched: " << today_iso() << " via geant4/g4_table_generator\n"
      << "# Columns: T_MeV, S_elec_MeV_cm2_per_g, S_nuc_MeV_cm2_per_g, "
      << "S_total_MeV_cm2_per_g, R_csda_g_per_cm2, R_proj_g_per_cm2, "
      << "detour_factor\n"
      << "T_MeV,S_elec_MeV_cm2_per_g,S_nuc_MeV_cm2_per_g,"
         "S_total_MeV_cm2_per_g,R_csda_g_per_cm2,R_proj_g_per_cm2,"
         "detour_factor\n";

  out << std::scientific << std::setprecision(6);
  for (double t_mev : ts_mev) {
    const double T = t_mev * MeV;
    // Stopping powers from G4EmCalculator: returns MeV/mm internally.
    // Convert by Geant4 unit system: divide by (MeV/(g/cm2)) i.e.
    //   S [MeV cm^2/g] = S_g4 / rho_g4  ; G4EmCalc has Compute*DEDX(rho-free)
    // Easier: use ComputeElectronicDEDX which returns MeV cm^2 / g when
    //   we divide by density. But G4EmCalculator returns G4 units, so
    //   convert via /(MeV/mm) and *( mm / cm) /(rho [g/cm3]).
    const double s_elec_g4 =
      calc.ComputeElectronicDEDX(T, part, mat);
    const double s_nuc_g4 =
      calc.ComputeNuclearDEDX(T, part, mat);
    // GetCSDARange requires G4EmParameters::SetBuildCSDARange(true)
    // *and* a particle/process that registers a CSDA table — not all
    // combinations build one. Fall back to GetRange (energy-loss range
    // table, always built) when the CSDA lookup returns zero.
    double r_csda_g4 = calc.GetCSDARange(T, part, mat);
    if (r_csda_g4 <= 0.0) {
      r_csda_g4 = calc.GetRange(T, part, mat);
    }

    // G4 stopping power is energy per length. Convert to MeV cm^2 / g
    // by dividing by density. Geant4 internal unit for energy/length is
    // (MeV)/(mm); to get MeV/cm multiply by 10. Then divide by rho.
    const double mev_per_cm_to_mass = 1.0 / rho_g_per_cm3;
    const double s_elec = s_elec_g4 / (MeV / mm) * 10.0 * mev_per_cm_to_mass;
    const double s_nuc = s_nuc_g4 / (MeV / mm) * 10.0 * mev_per_cm_to_mass;
    const double s_total = s_elec + s_nuc;
    // CSDA range in Geant4 internal units (mm). Grammage = R[cm] * rho.
    const double r_csda_cm = r_csda_g4 / mm * 0.1;
    const double r_csda = r_csda_cm * rho_g_per_cm3;

    out << t_mev << ',' << s_elec << ',' << s_nuc << ',' << s_total << ','
        << r_csda << ','
        // No separate projected range available; reuse CSDA. Detour=1.
        << r_csda << ',' << 1.0 << '\n';
  }

  out.close();
  std::cerr << "wrote " << args.output << " (" << ts_mev.size()
            << " rows) for " << args.particle << " in " << mat->GetName()
            << "\n";
  delete runManager;
  return 0;
}

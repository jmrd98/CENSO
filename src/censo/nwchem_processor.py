import os
import shutil
from collections import OrderedDict
from functools import reduce

from .utilities import od_insert, SolventHelper
from .logging import setup_logger
from .datastructure import GeometryData, ParallelJob
from .params import (
    CODING,
    USER_ASSETS_PATH,
    WARNLEN,
)
from .part import CensoPart

from .qm_processor import QmProc

logger = setup_logger(__name__)


class NWChemParser:
    """
    Parser for NWChem input and output files.
    Capable of reading input files and converting them to an OrderedDict,
    and writing input files from an OrderedDict.
    """

    def read_input(self, path: str) -> OrderedDict:
        """Reads NWChem input file and parses it into an OrderedDict."""
        with open(path, "r") as infile:
            lines = infile.readlines()
        return self.__todict(lines)

    def write_input(self, path: str, indict: OrderedDict) -> None:
        """Writes a NWChem input file from an OrderedDict."""
        with open(path, "w") as outfile:
            outfile.writelines(self.__tolines(indict))

    def __todict(self, lines: list[str]) -> OrderedDict:
        """Converts lines from NWChem input into an OrderedDict."""
        result = OrderedDict()
        current_section = None
        for line in lines:
            if line.strip().startswith("#"):
                continue  # Skip comments
            if line.strip().endswith(":"):
                current_section = line.strip().replace(":", "")
                result[current_section] = []
            elif current_section:
                result[current_section].append(line.strip())
        return result

    def __tolines(self, indict: OrderedDict) -> str:
        """
        Converts an OrderedDict back into a formatted NWChem input using Python's str.format() method.

        Args:
            indict (OrderedDict): Contains all sections like 'start', 'geometry', 'basis', etc.

        Returns:
            str: The complete NWChem input file content as a single string.
        """
        # Extract values from the OrderedDict
        job_name = indict.get("start", ["job"])[0]
        geometry = indict.get("geometry", [])[0]
        # print("geometry", geometry)
        # geometry is [{'element': 'O', 'xyz': [-2.0156358993, 0.4782335569, 0.2431094143]}, {'element': 'C', 'xyz': [-0.9402116723, 0.0302552922, -0.0526792532]}, {'element': 'C', 'xyz': [-0.725198202, -1.4427194194, -0.0897633946]}, {'element': 'H', 'xyz': [-1.0221728335, -1.8846780121, 0.8607375561]}, {'element': 'H', 'xyz': [-1.3304708138, -1.8929920856, -0.8755536896]}, {'element': 'H', 'xyz': [0.3220882644, -1.6883867431, -0.2705475036]}, {'element': 'C', 'xyz': [0.1878113278, 0.9550448642, -0.4060231877]}, {'element': 'H', 'xyz': [0.3196376861, 0.9600966566, -1.4917712337]}, {'element': 'H', 'xyz': [-0.099306823, 1.9728128902, -0.1231797244]}, {'element': 'C', 'xyz': [1.4678988806, 0.5660538741, 0.2408202311]}, {'element': 'H', 'xyz': [1.9474025289, 1.2315483811, 0.9510468963]}, {'element': 'O', 'xyz': [2.0142754, -0.4777864802, 0.0476423808]}]
        # needs to be in the format 'O 0.00000000 0.00000000 0.00000000\nC 0.75695031 -0.58587586 0.00000000\nH -0.75695031 -0.58587586 0.00000000\n'
        # with fixed field width
        geometry = "\n".join(
            [
                f"{atom['element']} {atom['xyz'][0]:.10f} {atom['xyz'][1]:.10f} {atom['xyz'][2]:.10f}"
                for atom in geometry
            ]
        )
        # print("geometry2", geometry)
        # exit(0)
        basis_set = "\n".join(indict.get("basis", []))
        dft_section = "\n".join(indict.get("dft", []))
        charge = indict.get("charge", ["0"])[0]  # Default charge is 0
        multiplicity = indict.get("multiplicity", ["1"])[0]  # Default multiplicity is 1
        task = indict.get("task", ["energy"])[0]
        solvation = "\n".join(indict.get("solvation", []))

        # Define the template using str.format() placeholders
        nwchem_template = """{job_name}
geometry units angstrom
{geometry} 
end
charge {charge}
{basis_set}
{solvation}
{dft_section}
{task}
"""

        # Format the template with the extracted values
        nwchem_input = nwchem_template.format(
            job_name=job_name,
            charge=charge,
            multiplicity=multiplicity,
            geometry=geometry.strip(),
            basis_set=basis_set,
            dft_section=dft_section,
            task=task,
            solvation=solvation,
        )
        # print("job_name", job_name)
        # print("charge", charge)
        # print("multiplicity", multiplicity)
        # print("geometry", geometry.strip())
        # print("basis_set", basis_set)
        # print("dft_section", dft_section)
        # print("task", task)
        # print("indict", indict)
        # print("nwchem_input", nwchem_input)
        # exit(0)
        return nwchem_input


class NWChemProc(QmProc):
    """
    NWChemProcessor class for managing NWChem calculations.
    Handles input creation, job execution, and output parsing for NWChem.
    """

    _progname = "nwchem"

    __cosmo_dcs = {
        "acetone": 20.7,
        "acetonitrile": 36.6,
        "aniline": 6.9,
        "benzaldehyde": 18.2,
        "benzene": 2.3,
        "ccl4": 2.2,
        "ch2cl2": 9.1,
        "chcl3": 4.8,
        "cs2": 2.6,
        "cyclohexane": 2.0,
        "dichloroethane": 10.125,
        "diethylether": 4.4,
        "dioxane": 2.2,
        "dmf": 38.3,
        "dmso": 47.2,
        "ethanol": 24.6,
        "ethylacetate": 5.9,
        "furan": 3.0,
        "h2o": 80.1,
        "hexadecane": 2.1,
        "hexane": 1.9,
        "methanol": 32.7,
        "nitromethane": 38.2,
        "octane": 1.94,
        "octanol": 9.9,
        "phenol": 8.0,
        "thf": 7.6,
        "toluene": 2.4,
        "woctanol": 8.1,
    }

    def __init__(self, *args, **kwargs):
        """
        Initialize the NWChemProcessor.

        Args:
            workdir (str): Directory where jobs are run.
            settings (dict): Dictionary containing settings for NWChem, such as basis set, charge, and NWChem path.
        """
        super().__init__(
            *args, **kwargs
        )  # Initialize QmProc with the working directory

        # Handle NWChem-specific initialization
        # self.settings = settings
        # self.nwchem_exec = settings.get(
        #     "nwchempath", "nwchem"
        # )  # Path to NWChem executable
        self._jobtypes = {
            **self._jobtypes,
            **{
                "sp": self._sp,
                "gsolv": self._gsolv,
                # "xtb_sp": self._xtb_sp,
                # "xtb_gsolv": self._xtb_gsolv,
                "opt": self._opt,
                "nmr": self._nmr,
                "uvvis": self._uvvis,
            },
        }
        # censopart = CensoPart()
        # self._get_settings = censopart.get_settings
        # self._get_general_settings = censopart.get_general_settings
        # Stores setting wether to copy MO-files for faster SCFs
        self.copy_mo: bool = False
        # Define the return code to error message mapping
        self.__returncode_to_err = {
            -1: "NWChem execution failed",
            1: "NWChem execution failed",
            2: "Input file error",
            3: "NWChem runtime error",
            4: "Convergence not achieved",
            # Add additional return codes and corresponding error messages as necessary
        }

    def _todict(self, lines: list[str]) -> OrderedDict:
        """
        Converts NWChem input lines into an OrderedDict.

        Args:
            lines: List of lines from the input file.

        Returns:
            OrderedDict: Structured representation of the input file.
        """
        result = OrderedDict()
        current_section = None
        for line in lines:
            if line.strip().startswith("#"):
                continue  # Skip comments
            if line.strip().endswith(":"):
                current_section = line.strip().replace(":", "")
                result[current_section] = []
            elif current_section:
                result[current_section].append(line.strip())
        return result

    def _tolines(self, indict: OrderedDict) -> list[str]:
        """
        Converts an OrderedDict back to NWChem input lines.

        Args:
            indict: OrderedDict representation of the input data.

        Returns:
            list[str]: List of lines to write to the input file.
        """
        lines = []
        for section, entries in indict.items():
            lines.append(f"{section}:\n")
            lines.extend([f"  {entry}\n" for entry in entries])
        return lines

    def _remove_comments(self, lines: list[str]) -> list[str]:
        """
        Removes comments from the input lines.

        Args:
            lines: List of lines from the input file.

        Returns:
            list[str]: Input lines with comments removed.
        """
        return [line for line in lines if not line.strip().startswith("#")]

    def _eob(self, indict: OrderedDict, key: str) -> bool:
        """
        Checks if the given key marks the end of a block.

        Args:
            indict: OrderedDict representation of the input data.
            key: Key to check.

        Returns:
            bool: True if the key marks the end of a block, False otherwise.
        """
        return key == "end"

    def _check_output(self, output_file: str) -> bool:
        """
        Check the NWChem output for errors or convergence issues.

        Args:
            output_file: Path to the NWChem output file.

        Returns:
            bool: True if the job completed successfully, False otherwise.
        """
        with open(output_file, "r") as f:
            for line in f:
                if "Total energy" in line:
                    return True
        return False

    def _apply_flags(self, job: ParallelJob, indict: OrderedDict) -> OrderedDict:
        """
        Apply specific flags to the NWChem input based on previous results.

        Args:
            job: ParallelJob object containing job information.
            indict: OrderedDict containing the current input data.

        Returns:
            OrderedDict: Updated input data with the applied flags.
        """
        if "scf_failed" in job.flags:
            indict["scf"] = ["maxiter 300", "thresh 1e-6"]
        return indict

    def _sp(
        self,
        job: ParallelJob,
        jobdir: str,
        filename="sp",
        no_solv: bool = False,
        prep: bool = True,
    ) -> tuple[dict[str, float | None], dict[str, any]]:
        """
        NWChem single-point calculation.

        Args:
            job: ParallelJob object containing the job information.
            jobdir: Path to the job directory.
            filename: Name of the input file.
            no_solv: If True, no solvent model is used.
            prep: If True, a new input file is generated.

        Returns:
            result: Dictionary containing the results of the calculation.
            meta: Metadata about the job.
        """
        # Initialize result and meta dictionaries
        result = {
            "energy": None,
        }
        meta = {
            "success": None,
            "error": None,
            "mo_path": None,
        }

        # Set input/output file paths
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # Prepare the input dictionary and write it to the input file
        if prep:
            indict = self.__prep(job, "sp", no_solv=no_solv)
            indict = self.__apply_flags(job, indict)

            # Write input file for NWChem
            parser = NWChemParser()
            parser.write_input(inputpath, indict)

        # Handle molecular orbitals if specified (if applicable to NWChem)
        # if self.copy_mo:
        #     self.__copy_mo(jobdir, filename, job.mo_guess)

        # Execute the NWChem job
        call = [f"{filename}.inp"]
        returncode, errors = self._make_call("nwchem", call, outputpath, jobdir)

        # Check if the job succeeded based on returncode
        meta["success"] = returncode == 0
        if not meta["success"]:
            logger.warning(f"Job for {job.name} failed. Stderr output:\n{errors}")

        # Parse the output to extract energy
        with open(outputpath, "r", encoding="utf-8") as out:
            lines = out.readlines()
        result["energy"] = next(
            (float(line.split()[4]) for line in lines if "Total DFT energy" in line),
            None,
        )

        # Check for errors in the output if the job succeeded
        if meta["success"]:
            meta["error"] = self.__check_output(lines)
            meta["success"] = meta["error"] is None and result["energy"] is not None
        else:
            meta["error"] = self.__returncode_to_err.get(returncode, "unknown_error")

        # Handle molecular orbitals if present (if applicable to NWChem)
        if self.copy_mo and os.path.isfile(os.path.join(jobdir, f"{filename}.movec")):
            meta["mo_path"] = os.path.join(jobdir, f"{filename}.movec")
        print("result", result)
        print("meta", meta)
        # exit()
        return result, meta

    def _opt(
        self,
        job: ParallelJob,
        jobdir: str,
        filename="opt",
        no_solv: bool = False,
        prep: bool = True,
    ) -> tuple[dict[str, float | None], dict[str, any]]:
        """
        NWChem geometry optimization.

        Args:
            job: ParallelJob object containing the job information.
            jobdir: Path to the job directory.
            filename: Name of the input file.
            no_solv: If True, no solvent model is used.
            prep: If True, a new input file is generated.

        Returns:
            result (dict[str, float | None]): Dictionary containing the results of the calculation.
            meta (dict[str, any]): Metadata about the job.
        """
        # Initialize result and meta dictionaries
        result = {
            "energy": None,
            "geometry": None,  # Store optimized geometry
        }
        meta = {
            "success": None,
            "error": None,
            "mo_path": None,
        }

        # Set input/output file paths
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # Prepare input for geometry optimization
        if prep:
            indict = self.__prep(job, "opt", no_solv=no_solv)
            indict = self.__apply_flags(job, indict)

            # Write input file for NWChem optimization
            parser = NWChemParser()
            parser.write_input(inputpath, indict)

        # Handle molecular orbitals if needed
        if self.copy_mo:
            self.__copy_mo(jobdir, filename, job.mo_guess)

        # Execute the NWChem job
        call = [f"{filename}.inp"]
        returncode, errors = self._make_call("nwchem", call, outputpath, jobdir)

        # Check if the job succeeded
        meta["success"] = returncode == 0
        if not meta["success"]:
            logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")

        # Parse the output for energy and geometry
        with open(outputpath, "r", encoding="utf-8") as out:
            lines = out.readlines()

        # Get final energy
        result["energy"] = next(
            (float(line.split()[4]) for line in lines if "Total SCF energy" in line),
            None,
        )

        # Extract optimized geometry
        geom_start = False
        geometry = []
        for line in lines:
            if "Output coordinates" in line:
                geom_start = True
                continue
            if geom_start and line.strip() == "":
                break
            if geom_start:
                geometry.append(line.strip())

        # Store the optimized geometry
        if geometry:
            result["geometry"] = "\n".join(geometry)

        # Check for errors and update metadata
        if meta["success"]:
            meta["error"] = self.__check_output(lines)
            meta["success"] = meta["error"] is None and result["energy"] is not None
        else:
            meta["error"] = self.__returncode_to_err.get(returncode, "unknown_error")

        # Handle molecular orbitals if present
        if self.copy_mo and os.path.isfile(os.path.join(jobdir, f"{filename}.movec")):
            meta["mo_path"] = os.path.join(jobdir, f"{filename}.movec")

        return result, meta

    # def _xtb_sp(
    #     self,
    #     job: ParallelJob,
    #     jobdir: str,
    #     filename: str = "xtb_sp",
    #     no_solv: bool = False,
    # ) -> tuple[dict[str, float | None], dict[str, any]]:
    #     """
    #     Calculates the single-point energy with GFNn-xTB or GFN-FF.

    #     Args:
    #         job (ParallelJob): job to run
    #         jobdir (str): path to the jobdir
    #         filename (str, optional): filename to use for the coord file. Defaults to "xtb_sp".
    #         no_solv (bool, optional): whether to run the sp in gas-phase. Defaults to False.

    #     Returns:
    #         result (dict[str, float | None]): result of the sp calculation
    #         meta (dict[str, any]): metadata about the job

    #     result = {
    #         "energy": None,
    #     }
    #     """
    #     # set results
    #     result = {
    #         "energy": None,
    #     }

    #     # set metadata
    #     meta = {
    #         "success": None,
    #         "error": None,
    #     }

    #     # set in/out path
    #     inputpath = os.path.join(jobdir, f"{filename}.coord")
    #     outputpath = os.path.join(jobdir, f"{filename}.out")
    #     xcontrolname = "xtb_sp-xcontrol-inp"

    #     # cleanup
    #     files = [
    #         "xtbrestart",
    #         "xtbtopo.mol",
    #         xcontrolname,
    #         "wbo",
    #         "charges",
    #         "gfnff_topo",
    #         f"{filename}.out",
    #     ]

    #     # remove potentially preexisting files to avoid confusion
    #     for file in files:
    #         if os.path.isfile(os.path.join(jobdir, file)):
    #             os.remove(os.path.join(jobdir, file))

    #     # generate coord file for xtb
    #     with open(inputpath, "w", newline=None) as file:
    #         file.writelines(job.conf.tocoord())

    #     # setup call for xtb single-point
    #     print("job", job)
    #     print(job.prepinfo)
    #     print()
    #     # exit()

    #     job.prepinfo["xtb_sp"] = {
    #         "gfnv": ,
    #         "solvent_key_xtb": SolventHelper.get_solvent(
    #             self._get_general_settings()["sm_rrho"],
    #             self._get_general_settings()["solvent"],
    #         ),
    #     }
    #     call = [
    #         f"{filename}.coord",
    #         "--" + job.prepinfo["xtb_sp"]["gfnv"],
    #         "--sp",
    #         "--chrg",
    #         f"{job.prepinfo['charge']}",
    #         "--norestart",
    #         "--parallel",
    #         f"{job.omp}",
    #     ]

    #     # add solvent to xtb call if not a gas-phase sp
    #     # (set either through run settings or by call kwarg e.g. for _xtb_gsolv)
    #     # NOTE on solvents_dict (or rather censo_solvents.json):
    #     # [0] is the normal name of the solvent, if it is available, [1] is the replacement
    #     if not (job.prepinfo["general"].get("gas-phase", False) or no_solv):
    #         call.extend(
    #             [
    #                 "--" + job.prepinfo["general"]["sm_rrho"],
    #                 job.prepinfo["xtb_sp"]["solvent_key_xtb"],
    #                 "reference",
    #                 "-I",
    #                 xcontrolname,
    #             ]
    #         )

    #         # set gbsa grid
    #         with open(os.path.join(jobdir, xcontrolname), "w", newline=None) as xcout:
    #             xcout.write("$gbsa\n")
    #             xcout.write("  gbsagrid=tight\n")
    #             xcout.write("$end\n")

    #     # call xtb
    #     returncode, errors = self._make_call("xtb", call, outputpath, jobdir)

    #     # if returncode != 0 then some error happened in xtb
    #     # TODO - returncodes
    #     if returncode != 0:
    #         meta["success"] = False
    #         meta["error"] = "unknown_error"
    #         logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
    #         return result, meta

    #     # read energy from outputfile
    #     with open(outputpath, "r", encoding=CODING, newline=None) as outputfile:
    #         for line in outputfile.readlines():
    #             if "| TOTAL ENERGY" in line:
    #                 result["energy"] = float(line.split()[3])
    #                 meta["success"] = True
    #             # TODO - important - what to do if calculation not converged?

    #     # FIXME - right now the case meta["success"] = None might appear if "TOTAL ENERGY" is not found in outputfile
    #     return result, meta

    # def _xtb_gsolv(
    #     self, job: ParallelJob, jobdir: str
    # ) -> tuple[dict[str, float | None], dict[str, any]]:
    #     """
    #     Calculate additive GBSA or ALPB solvation using GFNn-xTB or GFN-FF.

    #     Args:
    #         job (ParallelJob): job to run
    #         jobdir (str): path to the jobdir

    #     Returns:
    #         result (dict[str, any]): result of the gsolv calculation

    #     result = {
    #         "gsolv": None,
    #         "energy_xtb_gas": None,
    #         "energy_xtb_solv": None,
    #     }
    #     """
    #     # what is returned in the end
    #     result = {
    #         "gsolv": None,
    #         "energy_xtb_gas": None,
    #         "energy_xtb_solv": None,
    #     }

    #     meta = {
    #         "success": None,
    #         "error": None,
    #     }

    #     # run gas-phase GFN single-point
    #     spres, spmeta = self._xtb_sp(job, jobdir, filename="gas", no_solv=True)
    #     if spmeta["success"]:
    #         result["energy_xtb_gas"] = spres["energy"]
    #     else:
    #         meta["success"] = False
    #         meta["error"] = spmeta["error"]
    #         return result, meta

    #     # run single-point in solution:
    #     # ''reference'' corresponds to 1\;bar of ideal gas and 1\;mol/L of liquid
    #     #   solution at infinite dilution,
    #     spres, spmeta = self._xtb_sp(job, jobdir, filename="solv")
    #     if spmeta["success"]:
    #         result["energy_xtb_solv"] = spres["energy"]
    #     else:
    #         meta["success"] = False
    #         meta["error"] = spmeta["error"]
    #         return result, meta

    #     # only reached if both gas-phase and solvated sp succeeded
    #     result["gsolv"] = result["energy_xtb_solv"] - result["energy_xtb_gas"]
    #     meta["success"] = True

    #     return result, meta

    def _gsolv(
        self,
        job: ParallelJob,
        jobdir: str,
        filename="gsolv",
        no_solv: bool = False,
        prep: bool = True,
    ) -> tuple[dict[str, float | None], dict[str, any]]:
        """
        Run a solvent model calculation using NWChem (e.g., COSMO).

        Args:
            job: ParallelJob object containing the job information.
            jobdir: Path to the job directory.
            filename: Name of the input file.
            no_solv: If True, no solvent model is used.
            prep: If True, a new input file is generated.

        Returns:
            result: Dictionary containing the energy results.
            meta: Metadata about the job.
        """
        # Initialize result and meta dictionaries
        result = {
            "energy": None,
        }
        meta = {
            "success": None,
            "error": None,
        }

        # Set input/output file paths
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # Prepare input for solvation model calculation
        if prep:
            indict = self.__prep(job, "sp", no_solv=no_solv)
            indict = self.__apply_flags(job, indict)

            # Write input file for NWChem solvent calculation
            parser = NWChemParser()
            parser.write_input(inputpath, indict)

        # Execute the NWChem job
        call = [f"{filename}.inp"]
        returncode, errors = self._make_call("nwchem", call, outputpath, jobdir)

        # Check if the job succeeded
        meta["success"] = returncode == 0
        if not meta["success"]:
            logger.warning(
                f"Solvent model job for {job.conf.name} failed. Stderr output:\n{errors}"
            )

        # Parse the output to extract energy
        with open(outputpath, "r", encoding="utf-8") as out:
            lines = out.readlines()
        result["energy"] = next(
            (float(line.split()[4]) for line in lines if "Total DFT energy" in line),
            None,
        )

        # Check for errors and update metadata
        if meta["success"]:
            meta["error"] = self.__check_output(lines)
            meta["success"] = meta["error"] is None and result["energy"] is not None
        else:
            meta["error"] = self.__returncode_to_err.get(returncode, "unknown_error")

        return result, meta

    def _nmr(
        self,
        job: ParallelJob,
        jobdir: str,
        filename="nmr",
        prep: bool = True,
    ) -> tuple[dict[str, float | None], dict[str, any]]:
        """
        Run an NMR calculation using NWChem.

        Args:
            job: ParallelJob object containing the job information.
            jobdir: Path to the job directory.
            filename: Name of the input file.
            prep: If True, a new input file is generated.

        Returns:
            result: Dictionary containing the NMR shieldings.
            meta: Metadata about the job.
        """
        # Initialize result and meta dictionaries
        result = {
            "nmr_shieldings": None,
        }
        meta = {
            "success": None,
            "error": None,
        }

        # Set input/output file paths
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # Prepare input for NMR calculation
        if prep:
            indict = self.__prep(job, "nmr")
            indict = self.__apply_flags(job, indict)

            # Write input file for NWChem NMR calculation
            parser = NWChemParser()
            parser.write_input(inputpath, indict)

        # Execute the NWChem job
        call = [f"{filename}.inp"]
        returncode, errors = self._make_call("nwchem", call, outputpath, jobdir)

        # Check if the job succeeded
        meta["success"] = returncode == 0
        if not meta["success"]:
            logger.warning(
                f"NMR job for {job.conf.name} failed. Stderr output:\n{errors}"
            )

        # Parse the output to extract NMR shieldings
        with open(outputpath, "r", encoding="utf-8") as out:
            lines = out.readlines()

        # Extract NMR shieldings (if present)
        result["nmr_shieldings"] = []
        for line in lines:
            if "NMR shielding tensor" in line:
                # Parse NMR shieldings (this is highly dependent on the NWChem output format)
                result["nmr_shieldings"].append(line.strip())

        # Check for errors and update metadata
        if meta["success"]:
            meta["error"] = self.__check_output(lines)
            meta["success"] = meta["error"] is None and bool(result["nmr_shieldings"])
        else:
            meta["error"] = self.__returncode_to_err.get(returncode, "unknown_error")

        return result, meta

    def _uvvis(
        self,
        job: ParallelJob,
        jobdir: str,
        filename="uvvis",
        prep: bool = True,
    ) -> tuple[dict[str, any], dict[str, any]]:
        """
        Run a UV-Visible spectroscopy calculation using NWChem.

        Args:
            job: ParallelJob object containing the job information.
            jobdir: Path to the job directory.
            filename: Name of the input file.
            prep: If True, a new input file is generated.

        Returns:
            result: Dictionary containing the UV-Vis excitation data.
            meta: Metadata about the job.
        """
        result = {
            "excitations": [],
        }
        meta = {
            "success": None,
            "error": None,
        }

        # Set input/output file paths
        inputpath = os.path.join(jobdir, f"{filename}.inp")
        outputpath = os.path.join(jobdir, f"{filename}.out")

        # Prepare input for UV-Vis calculation
        if prep:
            indict = self.__prep(job, "uvvis")
            indict = self.__apply_flags(job, indict)

            # Write input file for NWChem UV-Vis calculation
            parser = NWChemParser()
            parser.write_input(inputpath, indict)

        # Execute the NWChem job
        call = [f"{filename}.inp"]
        returncode, errors = self._make_call("nwchem", call, outputpath, jobdir)

        # Check if the job succeeded
        meta["success"] = returncode == 0
        if not meta["success"]:
            logger.warning(
                f"UV-Vis job for {job.conf.name} failed. Stderr output:\n{errors}"
            )

        # Parse the output for excitation data
        with open(outputpath, "r", encoding="utf-8") as out:
            lines = out.readlines()

        # Extract excitation wavelengths and oscillator strengths
        for line in lines:
            if "Excitation energy" in line:
                # Parse excitation wavelengths and oscillator strengths (NWChem-specific format)
                result["excitations"].append(line.strip())

        # Check for errors and update metadata
        if meta["success"]:
            meta["error"] = self.__check_output(lines)
            meta["success"] = meta["error"] is None and bool(result["excitations"])
        else:
            meta["error"] = self.__returncode_to_err.get(returncode, "unknown_error")

        return result, meta

    # def _xtb_opt(
    #     self, job: ParallelJob, jobdir: str, filename: str = "xtb_opt"
    # ) -> tuple[dict[str, any], dict[str, any]]:
    #     """
    #     Geometry optimization using ANCOPT and ORCA gradients.
    #     Note that solvation is handled here always implicitly.

    #     Args:
    #         job: ParallelJob object containing the job information, metadata is stored in job.meta
    #         jobdir: path to the job directory
    #         filename: name of the input file

    #     Returns:
    #         result (dict[str, any]): dictionary containing the results of the calculation
    #         meta (dict[str, any]): metadata about the job

    #     Keywords required in prepinfo:
    #     - optcycles
    #     - hlow
    #     - optlevel
    #     - macrocycles
    #     - constraints

    #     result = {
    #         "energy": None,
    #         "cycles": None,
    #         "converged": None,
    #         "ecyc": None,
    #         "grad_norm": None,
    #         "geom": None,
    #     }
    #     """
    #     # NOTE: some "intuitivity problems":
    #     # the geometry of the conformer is written into a coord file and also into a xyz-file to be used by orca
    #     # xtb then outputs a file with the optimized geometry as 'xtbopt.coord', which is then read into the conformer
    #     # to update it's geometry

    #     # prepare result
    #     # 'ecyc' contains the energies for all cycles, 'cycles' stores the number of required cycles
    #     # 'gncyc' contains the gradient norms for all cycles
    #     # 'energy' contains the final energy of the optimization (converged or unconverged)
    #     # 'geom' stores the optimized geometry in GeometryData.xyz format
    #     result = {
    #         "energy": None,
    #         "cycles": None,
    #         "converged": None,
    #         "ecyc": None,
    #         "grad_norm": None,
    #         "gncyc": None,
    #         "geom": None,
    #     }

    #     meta = {
    #         "success": None,
    #         "error": None,
    #         "mo_path": None,
    #     }

    #     xcontrolname = "xtb_opt-xcontrol-inp"

    #     files = [
    #         "xtbrestart",
    #         "xtbtopo.mol",
    #         xcontrolname,
    #         "wbo",
    #         "charges",
    #         "gfnff_topo",
    #     ]

    #     # remove potentially preexisting files to avoid confusion
    #     for file in files:
    #         if os.path.isfile(os.path.join(jobdir, file)):
    #             os.remove(os.path.join(jobdir, file))

    #     # write conformer geometry to coord file
    #     with open(os.path.join(jobdir, f"{filename}.coord"), "w", newline=None) as file:
    #         file.writelines(job.conf.tocoord())

    #     # write xyz-file for orca
    #     with open(os.path.join(jobdir, f"{filename}.xyz"), "w", newline=None) as file:
    #         file.writelines(job.conf.toxyz())

    #     # set orca input path
    #     inputpath = os.path.join(jobdir, f"{filename}.inp")

    #     # prepare input dict
    #     parser = NWChemParser()
    #     indict = self.__prep(job, "xtb_opt", xyzfile=f"{filename}.xyz")

    #     # apply flags
    #     indict = self.__apply_flags(job, indict)

    #     # write orca input into file "xtb_opt.inp" in a subdir created for the
    #     # conformer
    #     parser.write_input(inputpath, indict)

    #     # append some additional lines to the coord file for ancopt
    #     with open(
    #         os.path.join(jobdir, f"{filename}.coord"), "a", newline=None
    #     ) as newcoord:
    #         newcoord.writelines(
    #             [
    #                 "$external\n",
    #                 f"   orca input file= {filename}.inp\n",
    #                 f"   orca bin= {self._paths['orcapath']}\n",
    #                 "$end\n",
    #             ]
    #         )

    #     # prepare configuration file for ancopt (xcontrol file)
    #     with open(os.path.join(jobdir, xcontrolname), "w", newline=None) as out:
    #         out.write("$opt \n")
    #         if job.prepinfo["xtb_opt"]["macrocycles"]:
    #             out.write(f"maxcycle={job.prepinfo['xtb_opt']['optcycles']} \n")
    #             out.write(f"microcycle={job.prepinfo['xtb_opt']['optcycles']} \n")

    #         out.writelines(
    #             [
    #                 "average conv=true \n",
    #                 f"hlow={job.prepinfo['xtb_opt']['hlow']} \n",
    #                 "s6=30.00 \n",
    #                 "engine=lbfgs\n",
    #                 "$external\n",
    #                 f"   orca input file= {filename}.inp\n",
    #                 f"   orca bin= {self._paths['orcapath']} \n",
    #             ]
    #         )

    #         # Import constraints
    #         if job.prepinfo["xtb_opt"]["constraints"] is not None:
    #             with open(job.prepinfo["xtb_opt"]["constraints"], "r") as f:
    #                 lines = f.readlines()

    #             out.writelines(lines)

    #         out.write("$end \n")

    #     # check, if there is an existing .gbw file and copy it if option
    #     # 'copy_mo' is true
    #     if self.copy_mo:
    #         self.__copy_mo(jobdir, filename, job.mo_guess)

    #     # prepare xtb call
    #     call = [
    #         f"{filename}.coord",  # name of the coord file generated above
    #         "--opt",
    #         job.prepinfo["xtb_opt"]["optlevel"],
    #         "--orca",
    #         "-I",
    #         xcontrolname,
    #     ]

    #     # set path to the ancopt output file
    #     outputpath = os.path.join(jobdir, f"{filename}.out")

    #     # call xtb
    #     returncode, errors = self._make_call("xtb", call, outputpath, jobdir)

    #     # check if optimization finished without errors
    #     # NOTE: right now, not converging scfs are not handled because returncodes need to be implemented first
    #     if returncode != 0:
    #         meta["success"] = False
    #         meta["error"] = "unknown_error"
    #         logger.warning(f"Job for {job.conf.name} failed. Stderr output:\n{errors}")
    #         # TODO - the xtb returncodes should be handled
    #         return result, meta

    #     # read output
    #     with open(outputpath, "r", encoding=CODING, newline=None) as file:
    #         lines = file.readlines()

    #     result["ecyc"] = []
    #     result["cycles"] = 0

    #     # Substrings indicating error in xtb
    #     error_ind = [
    #         "external code error",
    #         "|grad| > 500, something is totally wrong!",
    #         "abnormal termination of xtb",
    #     ]

    #     # Check if xtb terminated normally (if there are any error indicators
    #     # in the output)
    #     meta["success"] = (
    #         False
    #         if next((x for x in lines if any(y in x for y in error_ind)), None)
    #         is not None
    #         else True
    #     )
    #     if not meta["success"]:
    #         meta["error"] = "unknown_error"
    #         return result, meta

    #     # check convergence
    #     if (
    #         next((True for x in lines if "GEOMETRY OPTIMIZATION CONVERGED" in x), None)
    #         is True
    #     ):
    #         result["converged"] = True
    #     elif (
    #         next((True for x in lines if "FAILED TO CONVERGE GEOMETRY" in x), None)
    #         is True
    #     ):
    #         result["converged"] = False

    #     # Get the number of cycles
    #     if result["converged"] is not None:
    #         # tmp is one of the values from the dict defined below
    #         tmp = {
    #             True: ("GEOMETRY OPTIMIZATION CONVERGED", 5),
    #             False: ("FAILED TO CONVERGE GEOMETRY", 7),
    #         }
    #         tmp = tmp[result["converged"]]

    #         result["cycles"] = int(
    #             next(x for x in lines if tmp[0] in x).split()[tmp[1]]
    #         )

    #         # Get energies for each cycle
    #         result["ecyc"].extend(
    #             float(line.split("->")[-1])
    #             for line in filter(lambda x: "av. E: " in x, lines)
    #         )

    #         # Get all gradient norms for evaluation
    #         result["gncyc"] = [
    #             float(line.split()[3])
    #             for line in filter(lambda x: " gradient norm " in x, lines)
    #         ]

    #         # Get the last gradient norm
    #         result["grad_norm"] = result["gncyc"][-1]

    #         # store the final energy of the optimization in 'energy'
    #         result["energy"] = result["ecyc"][-1]
    #         meta["success"] = True

    #     if self.copy_mo:
    #         # store the path to the current .gbw file for this conformer
    #         meta["mo_path"] = os.path.join(jobdir, f"{filename}.gbw")

    #     # read out optimized geometry and update conformer geometry with this
    #     job.conf.fromcoord(os.path.join(jobdir, "xtbopt.coord"))
    #     result["geom"] = job.conf.xyz

    #     # TODO - this might be a case where it would be reasonable to raise an
    #     # exception
    #     try:
    #         assert result["converged"] is not None
    #     except AssertionError:
    #         meta["success"] = False
    #         meta["error"] = "unknown_error"

    #     return result, meta

    def __prep(
        self, job: ParallelJob, jobtype: str, no_solv: bool = False
    ) -> OrderedDict:
        """
        Prepare the input for a specific job type in NWChem.

        Args:
            job: ParallelJob object containing the job information.
            jobtype: The type of calculation (e.g., 'sp', 'opt', 'gsolv').
            no_solv: If True, no solvent model is used.

        Returns:
            OrderedDict: Prepared input data for the calculation.
        """
        indict = OrderedDict()
        print("prepinfo", job.prepinfo)
        # exit()
        indict = self.__prep_main(job.prepinfo, indict, jobtype)
        # Add start directive based on job or system name
        indict["start"] = [
            f"start {job.conf.name if hasattr(job.conf, 'name') else 'job'}"
        ]

        # Prepare the geometry section
        indict["geometry"] = [
            job.conf.xyz,  # Assuming job.conf.xyz holds properly formatted XYZ coordinates
        ]
        # print(indict["geometry"])
        # exit()
        # Prepare the basis set (use default from settings or allow customization)
        basis_set = indict["main"][1]

        indict["basis"] = ["basis", f"* library {basis_set}", "end"]

        # Prepare the DFT section with flexibility for different functionals
        functional = indict["main"][0]
        indict["dft"] = [
            "dft",
            f"  xc {functional}",
            "end",
        ]

        # # Add solvent model if applicable and not explicitly disabled
        if "sm" in job.prepinfo[jobtype] and not no_solv:
            print("sm check", job.prepinfo[jobtype]["sm"])
        if (
            not no_solv
            and not job.prepinfo["general"]["gas-phase"]
            and "sm" in job.prepinfo[jobtype]
        ):
            # indict["cosmo"] = ["cosmo", f"  solvent {self.settings['solvent']}", "end"]
            solv_helper = SolventHelper()
            solv = solv_helper.get_solvent(
                sm=job.prepinfo[jobtype]["sm"], name=job.prepinfo["general"]["solvent"]
            )
            dielect = self.__cosmo_dcs[solv]
            indict["solvation"] = ["cosmo", f"  dielec {dielect}", "end"]
            print("solvent", solv, "DIELECTRIC", self.__cosmo_dcs[solv])
            exit()

        # Set up task based on jobtype (e.g., single-point energy, optimization, etc.)
        if jobtype == "sp":
            indict["task"] = ["task dft energy"]
        elif jobtype == "opt":
            indict["task"] = ["task dft optimize"]
        elif jobtype == "gsolv" and not no_solv:
            indict["task"] = ["task dft energy"]
        elif jobtype == "nmr":
            indict["task"] = ["task dft nmr"]
        elif jobtype == "uvvis":
            indict["task"] = ["task dft excitation"]

        return indict

    # After building the input, write it to a file
    # input_dict = __prep(job, "sp")
    # write_nwchem_input("h2o.nw", input_dict)

    def __prep_main(
        self, prepinfo: dict[str, any], indict: OrderedDict, jobtype: str
    ) -> OrderedDict:
        if "main" not in indict:
            indict["main"] = []

        # grab func, basis
        # exit()
        if jobtype == "gsolv":
            print("jobtype", jobtype)
            jobtype = "xtb_gsolv"
            print("jobtype", jobtype)
            exit()
        # prepinfo["xtb_sp"] = {
        #     "gfnv": self.get_settings()["gfnv"],
        #     "solvent_key_xtb": SolventHelper.get_solvent(
        #         self.get_general_settings()["sm_rrho"],
        #         self.get_general_settings()["solvent"],
        #     ),
        # }
        print("prepinfo", prepinfo)
        print("jobtype", prepinfo[jobtype], jobtype)
        func = prepinfo[jobtype]["func_name"]
        basis = prepinfo[jobtype]["basis"]
        functype = prepinfo[jobtype]["func_type"]
        disp = prepinfo[jobtype]["disp"]

        indict["main"].append(func)
        indict["main"].append(basis)
        indict["main"].append(functype)
        indict["main"].append(disp)

        return indict

    def __apply_flags(self, job: ParallelJob, indict: OrderedDict) -> OrderedDict:
        """
        Apply specific flags to the NWChem input based on previous results.

        Args:
            job: ParallelJob object containing job information.
            indict: OrderedDict containing the current input data.

        Returns:
            OrderedDict: Updated input data with the applied flags.
        """
        if "scf_failed" in job.flags:
            indict["scf"] = ["maxiter 300", "thresh 1e-6"]

        return indict

    @staticmethod
    def __check_output(lines: list[str]) -> str | None:
        """
        Checks the lines from the output file for errors and returns them.

        Args:
            lines: list of lines from the output file.

        Returns:
            str | None: error message if an error was found, None otherwise
        """
        # Dict mapping specific messages from the output to error messages
        # TODO - this should be extended later
        out_to_err = {
            "SCF NOT CONVERGED": "scf_not_converged",
        }
        for line in lines:
            if any(key in line for key in out_to_err.keys()):
                # Returns the first error found
                key = next(filter(lambda x: x in line, out_to_err.keys()))
                return out_to_err[key]
        return None

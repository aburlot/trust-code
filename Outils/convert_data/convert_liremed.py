#!/usr/bin/env python3

"""
Update TRUST dataset to use new 'Lire_MED' syntax.

Authors: A Bruneton, E Saikali
"""
from tparser import TRUSTParser

class LireMEDConverter(TRUSTParser):
    """
    Former syntax: 
        Lire_MED [ vef|convertAllToPoly  ] [ family_names_from_group_names | short_family_names ] domaine_name mesh_name filename.med

    New syntax:
       lire_med {
          domaine dom
          file toto.med
          mesh the_mesh_in_file  // optional - if not there, first mesh taken.
          [convertAllToPoly]
          [exclure_groupes N GRP1 GRP2 ... GRPN ] // groups to exclude if any
       }
    """
    LIST_PARAMS = [("convertalltopoly", bool),
                  ("family_names_from_group_names", bool),
                  ("short_family_names", bool),
                  ("dom", str),
                  ("mesh", str),
                  ("file", str)
                   ]

    def __init__(self):
        super().__init__()
        self.lire_med = []

    def loadLireMED(self):
        ok, self.lire_med = self.loadNoCurlyBraceGeneric(["lire_med", "read_med"])
        return ok

    def createNewLireMED(self, it):
        dom, file, mesh, catp = it["dom"], it["file"], it["mesh"], it["convertalltopoly"]
        lm =  ["Lire_MED", "{", "\n"]
        lm += ["   domain", dom, "\n"]
        lm += ["   file", file, "\n"]
        if mesh != "--any--":
            lm += ["   mesh", mesh, "\n"]
        if catp:
            lm += ["   convertAllToPoly\n"]
        import os
        if not os.path.exists(file) :
            import warnings
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            warnings.warn(f"WARNING : {file} not found, exlude_groups in the new syntax of Read_MED block will be missing!!!!!")
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
        else :
            from trustutils import run
            run.useMEDCoupling()
            import medcoupling as mc
            a=mc.MEDFileData(file)
            m=a.getMeshes()[0]
            boundaries_from_grps = []

            for i in m.getFamiliesNames()[:] :
                nb_faces = m.getFamilyArr(-1,i,False).getNbOfElems()
                if nb_faces > 0 :
                    family_id = m.getFamilyId(i)
                    groups = m.getGroupsOnFamily(i)
                    if len(groups) == 1:
                        boundaries_from_grps.append(groups[0])
                    else:
                        for k in range(len(groups)):
                            nb_families = len(m.getFamiliesIdsOnGroup(groups[k]))
                            if nb_families == 1 :
                                boundaries_from_grps.append(groups[k])
            grpnames = m.getGroupsNames()
            grps_to_exclude = list(set(grpnames) - set(boundaries_from_grps))
            size = len(grps_to_exclude)
            grps_to_exclude=" ".join(grps_to_exclude)
            if size:
                lm += [f"   exclude_groups {size} {grps_to_exclude} \n"]

        lm += ["}\n"]
        return lm

    def outputData(self, fNameO):
        """ Write everything out.
        """
        tt = self.tabToken
        newData, prev = [], 0
        for lmb in self.lire_med:
            newData.extend(tt[prev:lmb["start"]])
            # New Lire_MED block
            lm = self.createNewLireMED(lmb)
            newData.extend(lm)
            prev = lmb["end"]
        # Finish writing:
        newData.extend(tt[prev:])
        self.unTokenizeAndWrite(newData, fNameO)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: convert_liremed.py <input_file.data> <output_file.data>")
        sys.exit(-1)
    fNameI, fNameO = sys.argv[1], sys.argv[2]
    dm = LireMEDConverter()
    dm.readAndTokenize(fNameI)
    if dm.loadLireMED():
        dm.outputData(fNameO)
        print("File '%s' written!" % fNameO)
    else:
        print("Nothing done.")
        sys.exit(-1)

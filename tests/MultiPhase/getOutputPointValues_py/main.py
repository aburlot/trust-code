# -*- coding: utf-8 -*-


class TRUSTICoCo():
    """
    TRUSTICoCoclass which implements the RunICoCo method
    """
    @classmethod
    def RunICoCo(self):
        """
        Execute ...
        """
        import medcoupling as mc
        import trusticoco as ti
        import time

        print ("Using ICoCo version",ti.ICOCO_VERSION)
        pbT = ti.ProblemTrio()
        pbT.name = "TRUST"
        pbT.setDataFile("jdd.data")
        pbT.initialize()

        porosity = mc.ICoCoMEDDoubleField()

        print("InputFieldsNames : ",pbT.getInputFieldsNames())
        print("OutputFieldsNames : ",pbT.getOutputFieldsNames())

        def run(pb):
            """
            Internal method for RunICoCo. It defines the time looping scheme.
            """
            stop = False # Does the Problem want to stop ?
            t = 0.0

            init = True
            while not stop:
                dt,stop = pbT.computeTimeStep()
                if stop: 
                    return

                pbT.initTimeStep(dt)
                
                x = ti.VecDouble()
                y = ti.VecDouble()
                z = ti.VecDouble()
                vals = ti.VecDouble()
                vals2 = ti.VecDouble()

                x.push_back(0.5)
                y.push_back(0.25)
                x.push_back(0.5)
                y.push_back(0.85)


                print("# --------------------------------------------------------------------")
                print("### Test Ref_champ VITESSE_LIQUIDE_SODIUM")
                pbT.getOutputPointValues("VITESSE_LIQUIDE_SODIUM", x,y,z, vals, 0)
                pbT.getOutputPointValues("VITESSE_LIQUIDE_SODIUM", x,y,z, vals2, 1)

                for i in range (2):
                    print ("@@@@@@@ VITESSE_LIQUIDE_SODIUM vals compo 0, vector : ", i, "  ", vals[i])

                for i in range (2):
                    print ("@@@@@@@ VITESSE_LIQUIDE_SODIUM vals compo 1, vector : ", i, "  ", vals2[i])

                a = pbT.getOutputPointValues("VITESSE_LIQUIDE_SODIUM", 0.5, 0.25, 0, 0)
                print ("@@@@@@@ VITESSE_LIQUIDE_SODIUM vals compo 0, point : ", a)
                a = pbT.getOutputPointValues("VITESSE_LIQUIDE_SODIUM", 0.5, 0.25, 0, 1)
                print ("@@@@@@@ VITESSE_LIQUIDE_SODIUM vals compo 1, point : ", a)

                print("# --------------------------------------------------------------------")
                print("### Test Ref_champ VITESSE_GAZ_SODIUM")
                pbT.getOutputPointValues("VITESSE_GAZ_SODIUM", x,y,z, vals, 0)
                pbT.getOutputPointValues("VITESSE_GAZ_SODIUM", x,y,z, vals2, 1)

                for i in range (2):
                    print ("@@@@@@@ VITESSE_GAZ_SODIUM vals compo 0, vector : ", i, "  ", vals[i])

                for i in range (2):
                    print ("@@@@@@@ VITESSE_GAZ_SODIUM vals compo 1, vector : ", i, "  ", vals2[i])

                a = pbT.getOutputPointValues("VITESSE_GAZ_SODIUM", 0.5, 0.25, 0, 0)
                print ("@@@@@@@ VITESSE_GAZ_SODIUM vals compo 0, point : ", a)
                a = pbT.getOutputPointValues("VITESSE_GAZ_SODIUM", 0.5, 0.25, 0, 1)
                print ("@@@@@@@ VITESSE_GAZ_SODIUM vals compo 1, point : ", a)
           
                print("# --------------------------------------------------------------------")
                print("### Test Ref_champ temperature")
                pbT.getOutputPointValues("temperature", x,y,z, vals, 0)
                pbT.getOutputPointValues("temperature", x,y,z, vals2, 1)

                for i in range (2):
                    print ("@@@@@@@ Temperature vals compo 0, vector : ", i, "  ", vals[i])

                for i in range (2):
                    print ("@@@@@@@ Temperature vals compo 1, vector : ", i, "  ", vals2[i])

                a = pbT.getOutputPointValues("temperature", 0.5, 0.25, 0, 0)
                print ("@@@@@@@ Temperature vals compo 0, point : ", a)
                a = pbT.getOutputPointValues("temperature", 0.5, 0.25, 0, 1)
                print ("@@@@@@@ Temperature vals compo 1, point : ", a)

                print("# --------------------------------------------------------------------")
                print("### Test Ref_champ TEMPERATURE_LIQUIDE_SODIUM & TEMPERATURE_GAZ_SODIUM")
                pbT.getOutputPointValues("TEMPERATURE_LIQUIDE_SODIUM", x,y,z, vals, 0)
                pbT.getOutputPointValues("TEMPERATURE_GAZ_SODIUM", x,y,z, vals2, 1)

                for i in range (2):
                    print ("@@@@@@@ temperature_liquide_sodium vals compo 0, vector : ", i, "  ", vals[i])

                for i in range (2):
                    print ("@@@@@@@ temperature_gaz_sodium vals compo 1, vector : ", i, "  ", vals2[i])

                a = pbT.getOutputPointValues("temperature_liquide_sodium", 0.5, 0.25, 0, 0)
                print ("@@@@@@@ temperature_liquide_sodium vals compo 0, point : ", a)
                a = pbT.getOutputPointValues("temperature_gaz_sodium", 0.5, 0.25, 0, 1)
                print ("@@@@@@@ temperature_gaz_sodium vals compo 1, point : ", a)

                print("# --------------------------------------------------------------------")
                print("### Test Ref_champ Alpha")
                pbT.getOutputPointValues("Alpha", x,y,z, vals, 0)
                pbT.getOutputPointValues("Alpha", x,y,z, vals2, 1)

                for i in range (2):
                    print ("@@@@@@@ Alpha vals compo 0, vector : ", i, "  ", vals[i])

                for i in range (2):
                    print ("@@@@@@@ Alpha vals compo 1, vector : ", i, "  ", vals2[i])

                a = pbT.getOutputPointValues("Alpha", 0.5, 0.25, 0, 0)
                print ("@@@@@@@ Alpha vals compo 0, point : ", a)
                a = pbT.getOutputPointValues("Alpha", 0.5, 0.25, 0, 1)
                print ("@@@@@@@ Alpha vals compo 1, point : ", a)

                print("# --------------------------------------------------------------------")
                print("### Test Champ_Generique H_ANALYTIQUE")
                pbT.getOutputPointValues("H_ANALYTIQUE", x,y,z, vals, 0)

                for i in range (2):
                    print ("@@@@@@@ H_ANALYTIQUE vals compo 0, vector : ", i, "  ", vals[i])

                a = pbT.getOutputPointValues("H_ANALYTIQUE", 0.5, 0.25, 0, 0)
                print ("@@@@@@@ H_ANALYTIQUE vals compo 0, point : ", a)


                print("# --------------------------------------------------------------------")
                print("### Test Champ_Generique ENTHALPIE_ELEM_dom")
                pbT.getOutputPointValues("ENTHALPIE_ELEM_dom", x,y,z, vals, 0)
                pbT.getOutputPointValues("ENTHALPIE_ELEM_dom", x,y,z, vals2, 1)

                for i in range (2):
                    print ("@@@@@@@ ENTHALPIE_ELEM_dom vals compo 0, vector : ", i, "  ", vals[i])

                for i in range (2):
                    print ("@@@@@@@ ENTHALPIE_ELEM_dom vals compo 1, vector : ", i, "  ", vals2[i])

                a = pbT.getOutputPointValues("ENTHALPIE_ELEM_dom", 0.5, 0.25, 0, 0)
                print ("@@@@@@@ ENTHALPIE_ELEM_dom vals compo 0, point : ", a)
                a = pbT.getOutputPointValues("ENTHALPIE_ELEM_dom", 0.5, 0.25, 0, 1)
                print ("@@@@@@@ ENTHALPIE_ELEM_dom vals compo 1, point : ", a)

                print("# --------------------------------------------------------------------")



                t = pbT.presentTime() # time av_totalanced

                # Main time loop:
                ok = pbT.solveTimeStep()
                t = pbT.presentTime() # time av_totalanced

                # iterate or stop ?
                if (not ok): # The resolution failed, try with a new time interval.
                    pbT.abortTimeStep()
                else: # The resolution was successful, validate and go to the next time step.
                    pbT.validateTimeStep()

                init = False
                print("=================================================================")

            stat = pbT.isStationary()
            if (stat):
                stop = True

        run(pbT)
        pbT.terminate()

if __name__ == "__main__":
    TRUSTICoCo().RunICoCo()

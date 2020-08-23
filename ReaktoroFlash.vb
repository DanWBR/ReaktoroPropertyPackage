'    Reaktoro Flash Algorithm
'    Copyright 2020 Daniel Wagner O. de Medeiros
'
'    This file is part of DWSIM.
'
'    DWSIM is free software: you can redistribute it and/or modify
'    it under the terms of the GNU General Public License as published by
'    the Free Software Foundation, either version 3 of the License, or
'    (at your option) any later version.
'
'    DWSIM is distributed in the hope that it will be useful,
'    but WITHOUT ANY WARRANTY; without even the implied warranty of
'    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
'    GNU General Public License for more details.
'
'    You should have received a copy of the GNU General Public License
'    along with DWSIM.  If not, see <http://www.gnu.org/licenses/>.

Imports DWSIM.Thermodynamics.PropertyPackages
Imports System.Math
Imports DWSIM.ExtensionMethods
Imports DWSIM.Thermodynamics
Imports DWSIM.Thermodynamics.PropertyPackages.Auxiliary.FlashAlgorithms
Imports DWSIM.Interfaces
Imports System.IO
Imports Python.Runtime
Imports DWSIM.Thermodynamics.BaseClasses
Imports DWSIM.GlobalSettings

<System.Serializable()> Public Class ReaktoroFlash

    Inherits FlashAlgorithm

    Public CompoundMaps As CompoundMapper

    Public Sub New()

        MyBase.New()

        If CompoundMaps Is Nothing Then CompoundMaps = New CompoundMapper()

    End Sub

    Public Overrides ReadOnly Property InternalUseOnly As Boolean
        Get
            Return True
        End Get
    End Property

    Public proppack As PropertyPackage

    Private P, Hf As Double

    Public Property CompoundProperties As List(Of ICompoundConstantProperties)

    Private Vxl0, Vf0 As Double()

    Private LoopVarF, LoopVarX As Double, LoopVarVz As Double(), LoopVarState As State

    Public Overloads Function Flash_PT(Vx As Array, T As Double, P As Double) As Dictionary(Of String, Object)

        Dim n As Integer = CompoundProperties.Count - 1
        Dim activcoeff(n) As Double
        Dim i As Integer

        'Vnf = feed molar amounts (considering 1 mol of feed)
        'Vnl = liquid phase molar amounts
        'Vnv = vapor phase molar amounts
        'Vns = solid phase molar amounts
        'Vxl = liquid phase molar fractions
        'Vxv = vapor phase molar fractions
        'Vxs = solid phase molar fractions
        'V, S, L = phase molar amounts (F = 1 = V + S + L)
        Dim K(n), Vnf(n), Vnl(n), Vnl_ant(n), Vxl(n), Vns(n), Vxs(n), Vnv(n), Vxv(n), Vxv_ant(n), Vf(n), V, S, L, Vp(n) As Double
        Dim sumN As Double = 0

        Dim saltonly As Boolean = False
        For i = 0 To Vx.Length - 1
            If Vx(i) > 0 And CompoundProperties(i).IsSalt Then
                saltonly = True
            ElseIf Vx(i) > 0 And Not CompoundProperties(i).IsSalt Then
                saltonly = False
                Exit For
            End If
        Next

        If saltonly Then

            'return flash calculation results.

            Dim resultsS As New Dictionary(Of String, Object)

            resultsS.Add("MixtureMoleFlows", Vx)
            resultsS.Add("VaporPhaseMoleFraction", 0.0)
            resultsS.Add("LiquidPhaseMoleFraction", 0.0)
            resultsS.Add("SolidPhaseMoleFraction", 1.0)
            resultsS.Add("VaporPhaseMolarComposition", proppack.RET_NullVector)
            resultsS.Add("LiquidPhaseMolarComposition", proppack.RET_NullVector)
            resultsS.Add("SolidPhaseMolarComposition", Vx)
            resultsS.Add("LiquidPhaseActivityCoefficients", proppack.RET_UnitaryVector)
            resultsS.Add("MoleSum", 1.0)

            Return resultsS

        End If

        Vnf = Vx.Clone

        Vf0 = Vx.Clone

        Dim names = proppack.RET_VNAMES().ToList
        Dim formulas As New List(Of String)

        For Each na In names
            If Not CompoundMaps.Maps.ContainsKey(na) Then
                Throw New Exception(String.Format("Compound {0} is not supported by this Property Package [{1}].", na, proppack.ComponentName))
            End If
        Next

        If Not PythonPathSet Then

            Dim ppath As String = Path.Combine(Path.GetDirectoryName(Reflection.Assembly.GetExecutingAssembly().Location), "reaktoro_python")
            Dim append As String = ppath + ";" + Path.Combine(ppath, "Library", "bin") + ";"

            Dim p1 As String = append + Environment.GetEnvironmentVariable("PATH", EnvironmentVariableTarget.Machine)
            ' Set Path
            Environment.SetEnvironmentVariable("PATH", p1, EnvironmentVariableTarget.Process)
            ' Set PythonHome
            Environment.SetEnvironmentVariable("PYTHONHOME", ppath, EnvironmentVariableTarget.Process)
            ' Set PythonPath
            Environment.SetEnvironmentVariable("PYTHONPATH", Path.Combine(p1, "Lib"), EnvironmentVariableTarget.Process)

            PythonPathSet = True

        End If

        If Not Settings.PythonInitialized Then

            proppack.Flowsheet.RunCodeOnUIThread(Sub()
                                                     PythonEngine.Initialize()
                                                     PythonEngine.BeginAllowThreads()
                                                 End Sub)

            Settings.PythonInitialized = True

        End If

        Dim speciesPhases As New Dictionary(Of String, String)
        Dim speciesAmounts As New Dictionary(Of String, Double)
        Dim speciesAmountsFinal As New Dictionary(Of String, Double)
        Dim compoundAmountsFinal As New Dictionary(Of String, Double)
        Dim inverseMaps As New Dictionary(Of String, String)

        Dim aqueous As String = "", gaseous As String = ""

        i = 0
        For Each na In names
            formulas.Add(CompoundMaps.Maps(na).Formula)
            speciesAmounts.Add(CompoundMaps.Maps(na).Formula, Vx(i))
            If CompoundMaps.Maps(na).AqueousName <> "" Then
                aqueous += CompoundMaps.Maps(na).AqueousName + " "
                speciesPhases.Add(CompoundMaps.Maps(na).AqueousName, "L")
                inverseMaps.Add(CompoundMaps.Maps(na).AqueousName, CompoundMaps.Maps(na).Formula)
            End If
            If CompoundMaps.Maps(na).GaseousName <> "" Then
                gaseous += CompoundMaps.Maps(na).GaseousName + " "
                speciesPhases.Add(CompoundMaps.Maps(na).GaseousName, "V")
                inverseMaps.Add(CompoundMaps.Maps(na).GaseousName, CompoundMaps.Maps(na).Formula)
            End If
            i += 1
        Next
        aqueous = aqueous.TrimEnd()
        gaseous = gaseous.TrimEnd()

        Using (Py.GIL())

            Dim reaktoro As Object = Py.Import("reaktoro")

            'Initialize a thermodynamic database
            Dim db = reaktoro.Database("supcrt98.xml")

            'Define the chemical system
            Dim editor = reaktoro.ChemicalEditor(db)

            editor.addAqueousPhase(aqueous)
            editor.addGaseousPhase(gaseous)

            'Construct the chemical system
            Dim mySystem = reaktoro.ChemicalSystem(editor)

            'Define the chemical equilibrium problem
            Dim problem = reaktoro.EquilibriumProblem(mySystem)
            problem.setTemperature(T, "kelvin")
            problem.setPressure(P, "pascal")

            For Each item In speciesAmounts
                problem.add(item.Key, item.Value, "mol")
            Next

            'Calculate the chemical equilibrium state
            Dim state = reaktoro.equilibrate(problem)

            Dim Ln As Double = state.phaseAmount("Aqueous").ToString().ToDoubleFromCurrent()
            Dim Vn As Double = state.phaseAmount("Gaseous").ToString().ToDoubleFromCurrent()
            Dim Sn As Double = 0.0

            Dim properties = state.properties

            Dim species = mySystem.species()
            Dim amounts = state.speciesAmounts()

            i = 0
            For Each item In species
                Dim name = item.name.ToString()
                speciesAmountsFinal.Add(name, amounts(i).ToString().ToDoubleFromCurrent())
                If Not compoundAmountsFinal.ContainsKey(inverseMaps(name)) Then
                    compoundAmountsFinal.Add(inverseMaps(name), 0.0)
                End If
                compoundAmountsFinal(inverseMaps(name)) += amounts(i).ToString().ToDoubleFromCurrent()
                If CompoundProperties(formulas.IndexOf(inverseMaps(name))).IsSalt Then
                    speciesPhases(name) = "S"
                    Sn += amounts(i).ToString().ToDoubleFromCurrent()
                    Ln -= amounts(i).ToString().ToDoubleFromCurrent()
                End If
                i += 1
            Next

            i = 0
            For Each item In species
                Dim name = item.name.ToString()
                Dim index = formulas.IndexOf(inverseMaps(name))
                Select Case speciesPhases(name)
                    Case "V"
                        Vxv(index) = amounts(i).ToString().ToDoubleFromCurrent()
                    Case "L"
                        Vxl(index) = amounts(i).ToString().ToDoubleFromCurrent()
                    Case "S"
                        Vxs(index) = amounts(i).ToString().ToDoubleFromCurrent()
                End Select
                Vnf(index) = compoundAmountsFinal(inverseMaps(name))
                i += 1
            Next

            Vxv = Vxv.NormalizeY()
            Vxl = Vxl.NormalizeY()
            Vxs = Vxs.NormalizeY()

            Dim ac = properties.lnActivityCoefficients().val

            i = 0
            For Each item In ac
                Dim index As Integer = formulas.IndexOf(inverseMaps(species(i).name.ToString()))
                activcoeff(index) = Math.Exp(item.ToString().ToDoubleFromCurrent())
                i += 1
            Next

            V = Vn / (Vn + Ln + Sn)
            L = Ln / (Vn + Ln + Sn)
            S = Sn / (Vn + Ln + Sn)

            sumN = Vn + Ln + Sn

        End Using

        'return flash calculation results.

        Dim results As New Dictionary(Of String, Object)

        results.Add("MixtureMoleFlows", Vnf)
        results.Add("VaporPhaseMoleFraction", V)
        results.Add("LiquidPhaseMoleFraction", L)
        results.Add("SolidPhaseMoleFraction", S)
        results.Add("VaporPhaseMolarComposition", Vxv)
        results.Add("LiquidPhaseMolarComposition", Vxl)
        results.Add("SolidPhaseMolarComposition", Vxs)
        results.Add("LiquidPhaseActivityCoefficients", activcoeff)
        results.Add("MoleSum", sumN)

        Return results

    End Function

    Public Overloads Function Flash_PH(ByVal Vz As Double(), ByVal P As Double, ByVal H As Double, ByVal Tref As Double) As Dictionary(Of String, Object)

        Dim Vn(1) As String, Vx(1), Vy(1), Vx_ant(1), Vy_ant(1), Vp(1), Ki(1), Ki_ant(1), fi(1) As Double
        Dim n, ecount As Integer
        Dim d1, d2 As Date, dt As TimeSpan
        Dim L, V, T, Pf As Double

        d1 = Date.Now

        n = Vz.Length - 1

        Vf0 = Vz.Clone

        Hf = H
        Pf = P

        ReDim Vn(n), Vx(n), Vy(n), Vx_ant(n), Vy_ant(n), Vp(n), Ki(n)

        If Tref = 0.0# Then Tref = 298.15
        If Double.IsNaN(Tref) Then Tref = 298.15

        Me.P = P

        Dim cnt As Integer = 0

        Dim fx, fx2, dfdx, x1, x0, dx As Double

        x1 = Tref

        Do

            If cnt < 2 Then

                fx = Herror({x1})
                fx2 = Herror({x1 + 0.1})

                dfdx = (fx2 - fx) / 0.1

            Else

                fx2 = fx
                fx = Herror({x1})

                dfdx = (fx - fx2) / (x1 - x0)

            End If

            If Abs(fx) <= 0.01 Then Exit Do

            dx = fx / dfdx

            x0 = x1
            x1 = x1 - dx

            If Double.IsNaN(x1) Or cnt > 25 Then
                Throw New Exception("PH Flash [Electrolyte]: Invalid result: Temperature did not converge.")
            End If

            cnt += 1

        Loop

        T = x1

        Dim tmp As Object = Flash_PT(Vz, T, P)

        Dim S, Vs(), Vnf(), sumN As Double

        sumN = tmp("MoleSum")
        L = tmp("LiquidPhaseMoleFraction")
        V = tmp("VaporPhaseMoleFraction")
        S = tmp("SolidPhaseMoleFraction")
        Vx = tmp("LiquidPhaseMolarComposition")
        Vy = tmp("VaporPhaseMolarComposition")
        Vs = tmp("SolidPhaseMolarComposition")
        Vnf = tmp("MixtureMoleFlows")

        d2 = Date.Now

        dt = d2 - d1

        WriteDebugInfo("PH Flash [Electrolyte]: Converged in " & ecount & " iterations. Time taken: " & dt.TotalMilliseconds & " ms.")

        'return flash calculation results.

        Dim results As New Dictionary(Of String, Object)

        results.Add("MixtureMoleFlows", Vnf)
        results.Add("VaporPhaseMoleFraction", V)
        results.Add("LiquidPhaseMoleFraction", L)
        results.Add("SolidPhaseMoleFraction", S)
        results.Add("VaporPhaseMolarComposition", Vy)
        results.Add("LiquidPhaseMolarComposition", Vx)
        results.Add("SolidPhaseMolarComposition", Vs)
        results.Add("MoleSum", sumN)
        results.Add("Temperature", T)
        results.Add("LiquidPhaseActivityCoefficients", tmp("LiquidPhaseActivityCoefficients"))

        Return results

    End Function

    Function Herror(ByVal x() As Double) As Double
        Return OBJ_FUNC_PH_FLASH(x(0), P, Vf0.Clone)
    End Function

    Function OBJ_FUNC_PH_FLASH(ByVal T As Double, ByVal P As Double, ByVal Vz As Object) As Double

        Dim tmp As Dictionary(Of String, Object) = Flash_PT(Vz, T, P)

        Dim FW0, FW, L, V, S, Vx(), Vy(), Vs(), sumN, _Hv, _Hl, _Hs As Double

        Dim n = Vz.Length - 1

        sumN = tmp("MoleSum")
        L = tmp("LiquidPhaseMoleFraction")
        V = tmp("VaporPhaseMoleFraction")
        S = tmp("SolidPhaseMoleFraction")
        Vx = tmp("LiquidPhaseMolarComposition")
        Vy = tmp("VaporPhaseMolarComposition")
        Vs = tmp("SolidPhaseMolarComposition")

        _Hv = 0
        _Hl = 0
        _Hs = 0

        Dim mmm, mmg, mml, mms As Double

        If V > 0.0# Then _Hv = proppack.DW_CalcEnthalpy(Vy, T, P, State.Vapor)
        If L > 0.0# Then _Hl = proppack.DW_CalcEnthalpy(Vx, T, P, State.Liquid)
        If S > 0.0# Then _Hs = proppack.DW_CalcSolidEnthalpy(T, Vs, CompoundProperties)

        mmg = proppack.AUX_MMM(Vy)
        mml = proppack.AUX_MMM(Vx)
        mms = proppack.AUX_MMM(Vs)

        mmm = V * mmg + L * mml + S * mml

        FW0 = 0.001 * proppack.AUX_MMM(Vz) 'kg
        FW = 0.001 * sumN * mmm 'kg

        Dim herr As Double = FW0 * Hf - FW * (((mmg * V / (mmg * V + mml * L + mms * S)) * _Hv + (mml * L / (mmg * V + mml * L + mms * S)) * _Hl + (mms * S / (mmg * V + mml * L + mms * S)) * _Hs))

        Return herr

        WriteDebugInfo("PH Flash [Electrolyte]: Current T = " & T & ", Current H Error = " & herr)

    End Function

    Private Function EnthalpyTx(ByVal x As Double, ByVal otherargs As Object) As Double

        Dim er As Double = LoopVarF - proppack.DW_CalcEnthalpy(LoopVarVz, x, LoopVarX, LoopVarState)
        Return er

    End Function

    Public Overloads Function Flash_TV(ByVal Vz As Double(), ByVal T As Double, ByVal V As Double, ByVal Pref As Double) As Object

        Dim n, ecount As Integer
        Dim d1, d2 As Date, dt As TimeSpan

        d1 = Date.Now

        Dim maxitINT As Integer = 100
        Dim maxitEXT As Integer = 100
        Dim tolINT As Double = 0.00001
        Dim tolEXT As Double = 0.00001

        n = Vz.Length - 1

        Dim Vx(n), Vy(n), Vp(n), Vcalc, Vspec, P, x, x0, x00, fx, fx0, fx00, Pmin, Pmax As Double

        Dim nl As New NestedLoops
        Dim flashresult = nl.CalculateEquilibrium(FlashSpec.T, FlashSpec.VAP, T, 0.0#, proppack, Vz, Nothing, Pref)
        Pmax = flashresult.CalculatedPressure
        flashresult = nl.CalculateEquilibrium(FlashSpec.T, FlashSpec.VAP, T, 1.0#, proppack, Vz, Nothing, Pref)
        Pmin = flashresult.CalculatedPressure

        P = Pmin + (1 - V) * (Pmax - Pmin)

        ecount = 0
        Vspec = V
        x = P

        Dim tmp As Dictionary(Of String, Object)

        Do

            tmp = Flash_PT(Vz, T, x)

            Vcalc = tmp("VaporPhaseMoleFraction")

            fx00 = fx0
            fx0 = fx

            fx = Vspec - Vcalc

            If Abs(fx) < tolEXT Then Exit Do

            x00 = x0
            x0 = x

            If ecount <= 1 Then
                x *= 0.99
            Else
                x = x - fx * (x - x00) / (fx - fx00)
                If Double.IsNaN(x) Then Throw New Exception(Calculator.GetLocalString("PropPack_FlashError"))
            End If

            ecount += 1

            If ecount > maxitEXT Then Throw New Exception(Calculator.GetLocalString("PropPack_FlashMaxIt2"))

        Loop

        P = x

        d2 = Date.Now

        dt = d2 - d1

        WriteDebugInfo("TV Flash [Sour Water]: Converged in " & ecount & " iterations. Time taken: " & dt.TotalMilliseconds & " ms.")

        Dim results As New Dictionary(Of String, Object)

        results.Add("MixtureMoleFlows", tmp("MixtureMoleFlows"))
        results.Add("VaporPhaseMoleFraction", tmp("VaporPhaseMoleFraction"))
        results.Add("LiquidPhaseMoleFraction", tmp("LiquidPhaseMoleFraction"))
        results.Add("SolidPhaseMoleFraction", tmp("SolidPhaseMoleFraction"))
        results.Add("VaporPhaseMolarComposition", tmp("VaporPhaseMolarComposition"))
        results.Add("LiquidPhaseMolarComposition", tmp("LiquidPhaseMolarComposition"))
        results.Add("SolidPhaseMolarComposition", tmp("SolidPhaseMolarComposition"))
        results.Add("MoleSum", tmp("MoleSum"))
        results.Add("Pressure", P)
        results.Add("LiquidPhaseActivityCoefficients", tmp("LiquidPhaseActivityCoefficients"))

        Return results

    End Function

    Public Overloads Function Flash_PV(ByVal Vz As Double(), ByVal P As Double, ByVal V As Double, ByVal Tref As Double) As Object

        Dim n, ecount As Integer
        Dim d1, d2 As Date, dt As TimeSpan

        d1 = Date.Now

        Dim maxitINT As Integer = 100
        Dim maxitEXT As Integer = 100
        Dim tolINT As Double = 0.00001
        Dim tolEXT As Double = 0.00001

        n = Vz.Length - 1

        Dim Vx(n), Vy(n), Vp(n), gamma(n), Vcalc, Vspec, T, x, x0, x00, fx, fx0, fx00, Pcalc As Double

        Dim nl As New NestedLoops

        ecount = 0
        Vspec = V
        x = Tref

        Dim tmp As Dictionary(Of String, Object)

        If Vspec = 0.0# Then

            'bubble point

            Do

                Vp = proppack.RET_VPVAP(x)

                tmp = Flash_PT(Vz, x, P)

                gamma = tmp("LiquidPhaseActivityCoefficients")

                fx00 = fx0
                fx0 = fx

                Pcalc = 0.0#
                For i = 0 To n
                    If Not CompoundProperties(i).IsIon And Not CompoundProperties(i).IsSalt Then
                        Pcalc += gamma(i) * Vp(i)
                    End If
                Next

                fx = P - Pcalc

                If Abs(fx) < 1.0# Then Exit Do

                x00 = x0
                x0 = x

                If ecount <= 1 Then
                    x += 1.0#
                Else
                    If ecount < 10 Then
                        x = x - 0.1 * fx * (x - x00) / (fx - fx00)
                    Else
                        x = x - fx * (x - x00) / (fx - fx00)
                    End If
                    If Double.IsNaN(x) Then Throw New Exception(Calculator.GetLocalString("PropPack_FlashError"))
                End If

                ecount += 1

                If ecount > 1000 Then Throw New Exception(Calculator.GetLocalString("PropPack_FlashMaxIt2"))

            Loop

        Else

            Do

                tmp = Flash_PT(Vz, x, P)

                Vcalc = tmp("VaporPhaseMoleFraction")

                fx00 = fx0
                fx0 = fx

                fx = Vspec - Vcalc

                If Abs(fx) < tolEXT Then Exit Do

                x00 = x0
                x0 = x

                If ecount <= 1 Then
                    x += 1.0#
                Else
                    x = x - 0.3 * fx * (x - x00) / (fx - fx00)
                    If Double.IsNaN(x) Then Throw New Exception(Calculator.GetLocalString("PropPack_FlashError"))
                End If

                ecount += 1

                If ecount > maxitEXT Then Throw New Exception(Calculator.GetLocalString("PropPack_FlashMaxIt2"))

            Loop

        End If

        T = x

        d2 = Date.Now

        dt = d2 - d1

        WriteDebugInfo("PV Flash [Sour Water]: Converged in " & ecount & " iterations. Time taken: " & dt.TotalMilliseconds & " ms.")

        Dim results As New Dictionary(Of String, Object)

        results.Add("MixtureMoleFlows", tmp("MixtureMoleFlows"))
        results.Add("VaporPhaseMoleFraction", tmp("VaporPhaseMoleFraction"))
        results.Add("LiquidPhaseMoleFraction", tmp("LiquidPhaseMoleFraction"))
        results.Add("SolidPhaseMoleFraction", tmp("SolidPhaseMoleFraction"))
        results.Add("VaporPhaseMolarComposition", tmp("VaporPhaseMolarComposition"))
        results.Add("LiquidPhaseMolarComposition", tmp("LiquidPhaseMolarComposition"))
        results.Add("SolidPhaseMolarComposition", tmp("SolidPhaseMolarComposition"))
        results.Add("MoleSum", tmp("MoleSum"))
        results.Add("Temperature", T)
        results.Add("LiquidPhaseActivityCoefficients", tmp("LiquidPhaseActivityCoefficients"))

        Return results

    End Function

    Public Function eval_h(ByVal n As Integer, ByVal x As Double(), ByVal new_x As Boolean, ByVal obj_factor As Double, ByVal m As Integer, ByVal lambda As Double(),
     ByVal new_lambda As Boolean, ByVal nele_hess As Integer, ByRef iRow As Integer(), ByRef jCol As Integer(), ByRef values As Double()) As Boolean

        If values Is Nothing Then

            Dim row(nele_hess - 1), col(nele_hess - 1) As Integer

            iRow = row
            jCol = col

        Else

        End If

        Return True

    End Function

    Public Overrides ReadOnly Property AlgoType As Enums.FlashMethod
        Get
            Return Enums.FlashMethod.Electrolyte
        End Get
    End Property

    Public Overrides ReadOnly Property Description As String
        Get
            Return "Electrolyte Flash"
        End Get
    End Property

    Public Overrides Function Flash_PH(Vz() As Double, P As Double, H As Double, Tref As Double, PP As PropertyPackage, Optional ReuseKI As Boolean = False, Optional PrevKi() As Double = Nothing) As Object
        proppack = PP
        Dim results = Flash_PH(Vz, P, H, Tref)
        'Return New Object() {L, V, Vx, Vy, T, ecount, Ki, 0.0#, PP.RET_NullVector, 0.0#, PP.RET_NullVector}
        With results
            Return New Object() {results("LiquidPhaseMoleFraction"), results("VaporPhaseMoleFraction"), results("LiquidPhaseMolarComposition"),
                                 results("VaporPhaseMolarComposition"), results("Temperature"), 1, PP.RET_NullVector, 0.0#, PP.RET_NullVector, results("SolidPhaseMoleFraction"), results("SolidPhaseMolarComposition")}
        End With
    End Function

    Public Overrides Function Flash_PS(Vz() As Double, P As Double, S As Double, Tref As Double, PP As PropertyPackage, Optional ReuseKI As Boolean = False, Optional PrevKi() As Double = Nothing) As Object
        Throw New NotImplementedException("Pressure-Entropy Flash Not Implemented")
    End Function

    Public Overrides Function Flash_PT(Vz() As Double, P As Double, T As Double, PP As PropertyPackage, Optional ReuseKI As Boolean = False, Optional PrevKi() As Double = Nothing) As Object
        proppack = PP
        Dim results = Flash_PT(Vz, T, P)
        With results
            Return New Object() {results("LiquidPhaseMoleFraction"), results("VaporPhaseMoleFraction"), results("LiquidPhaseMolarComposition"),
                                 results("VaporPhaseMolarComposition"), 1, 0.0#, PP.RET_NullVector, results("SolidPhaseMoleFraction"), results("SolidPhaseMolarComposition")}
        End With
    End Function

    Public Overrides Function Flash_PV(Vz() As Double, P As Double, V As Double, Tref As Double, PP As PropertyPackage, Optional ReuseKI As Boolean = False, Optional PrevKi() As Double = Nothing) As Object
        proppack = PP
        Dim results = Flash_PV(Vz, P, V, Tref)
        With results
            Return New Object() {results("LiquidPhaseMoleFraction"), results("VaporPhaseMoleFraction"), results("LiquidPhaseMolarComposition"),
                                 results("VaporPhaseMolarComposition"), results("Temperature"), 1, PP.RET_NullVector, 0.0#, PP.RET_NullVector, results("SolidPhaseMoleFraction"), results("SolidPhaseMolarComposition")}
        End With
    End Function

    Public Overrides Function Flash_TV(Vz() As Double, T As Double, V As Double, Pref As Double, PP As PropertyPackage, Optional ReuseKI As Boolean = False, Optional PrevKi() As Double = Nothing) As Object
        proppack = PP
        Dim results = Flash_TV(Vz, T, V, Pref)
        With results
            Return New Object() {results("LiquidPhaseMoleFraction"), results("VaporPhaseMoleFraction"), results("LiquidPhaseMolarComposition"),
                                 results("VaporPhaseMolarComposition"), results("Pressure"), 1, PP.RET_NullVector, 0.0#, PP.RET_NullVector, results("SolidPhaseMoleFraction"), results("SolidPhaseMolarComposition")}
        End With
    End Function

    Public Overrides ReadOnly Property Name As String
        Get
            Return "Electrolyte SVLE"
        End Get
    End Property

    Public Overrides ReadOnly Property MobileCompatible As Boolean
        Get
            Return False
        End Get
    End Property

End Class
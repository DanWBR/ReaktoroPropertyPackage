Imports System.IO
Imports DWSIM.Thermodynamics.PropertyPackages
Imports Python.Runtime
Imports DWSIM.ExtensionMethods
Imports DWSIM.GlobalSettings

Public Class ActivityCoefficients

    Public Overloads Function Calculate(Vx As Double(), T As Double, P As Double, pp As PropertyPackage) As Double()

        If Vx.SumY = 0.0 Then Return pp.RET_UnitaryVector()

        Dim i As Integer

        Dim CompProps = pp.DW_GetConstantProperties()
        Dim saltonly As Boolean = False
        For i = 0 To Vx.Length - 1
            If Vx(i) > 0 And CompProps(i).IsSalt Then
                saltonly = True
            ElseIf Vx(i) > 0 And Not CompProps(i).IsSalt Then
                saltonly = False
                Exit For
            End If
        Next

        If saltonly Then Return pp.RET_UnitaryVector()

        Dim CompoundMaps = New CompoundMapper()

        Dim n As Integer = Vx.Length - 1
        Dim activcoeff(n) As Double

        Dim names = pp.RET_VNAMES().ToList
        Dim formulas As New List(Of String)

        For Each na In names
            If Not CompoundMaps.Maps.ContainsKey(na) Then
                'Throw New Exception(String.Format("Compound {0} is not supported by this Property Package [{1}].", na, pp.ComponentName))
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

            AddDllDirectory(ppath)
            AddDllDirectory(Path.Combine(ppath, "Library", "bin"))

        End If

        If Not Settings.PythonInitialized Then

            pp.Flowsheet.RunCodeOnUIThread(Sub()
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

        Dim aqueous As String = ""

        i = 0
        For Each na In names
            formulas.Add(CompoundMaps.Maps(na).Formula)
            speciesAmounts.Add(CompoundMaps.Maps(na).Formula, Vx(i))
            If CompoundMaps.Maps(na).AqueousName <> "" Then
                aqueous += CompoundMaps.Maps(na).AqueousName + " "
                speciesPhases.Add(CompoundMaps.Maps(na).AqueousName, "L")
                inverseMaps.Add(CompoundMaps.Maps(na).AqueousName, CompoundMaps.Maps(na).Formula)
            End If
            i += 1
        Next
        aqueous = aqueous.TrimEnd()

        Dim pystate = Py.GIL()

        Dim ex0 As Exception = Nothing

        Try

            Dim sys As Object = PythonEngine.ImportModule("sys")

            Dim codeToRedirectOutput As String = "import sys" & Environment.NewLine + "from io import BytesIO as StringIO" & Environment.NewLine + "sys.stdout = mystdout = StringIO()" & Environment.NewLine + "sys.stdout.flush()" & Environment.NewLine + "sys.stderr = mystderr = StringIO()" & Environment.NewLine + "sys.stderr.flush()"

            PythonEngine.RunSimpleString(codeToRedirectOutput)

            Dim reaktoro As Object = Py.Import("reaktoro")
            Dim np As Object = Py.Import("numpy")

            'Initialize a thermodynamic database
            Dim db = reaktoro.Database("supcrt07-organics.xml")

            'Define the chemical system
            Dim editor = reaktoro.ChemicalEditor(db)

            editor.addAqueousPhase(aqueous)

            'Construct the chemical system
            Dim mySystem = reaktoro.ChemicalSystem(editor)

            Dim mols = np.fromiter(speciesAmounts.Values.ToArray(), np.float64)

            Dim props = reaktoro.ChemicalProperties(mySystem)
            props.update(T, P, mols)

            Dim species = mySystem.species()

            Dim ac = props.lnActivityCoefficients().val

            i = 0
            For Each item In ac
                Dim index As Integer = formulas.IndexOf(inverseMaps(species(i).name.ToString()))
                activcoeff(index) = item.ToString().ToDoubleFromInvariant()
                i += 1
            Next

        Catch ex As Exception

            pp.Flowsheet?.ShowMessage("Reaktoro error: " + ex.Message, DWSIM.Interfaces.IFlowsheet.MessageType.GeneralError)
            ex0 = ex

        Finally

            pystate?.Dispose()
            pystate = Nothing

        End Try

        If ex0 IsNot Nothing Then
            Throw ex0
        Else
            Return activcoeff.ExpY()
        End If

    End Function


End Class

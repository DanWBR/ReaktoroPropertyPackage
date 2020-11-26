﻿Imports FileHelpers
Imports System.IO

<DelimitedRecord(vbTab)> <IgnoreFirst()> <System.Serializable()>
Public Class SetschenowCoefficient

    Public Property Name As String
    <FieldNullValue(0.0)> Public Property Coefficient As String

End Class

Public Class SetschenowCoefficients

    Private Property Maps As New Dictionary(Of String, SetschenowCoefficient)

    Public Sub New()

        Dim map() As SetschenowCoefficient
        Dim fh1 As New FileHelperEngine(Of SetschenowCoefficient)

        Dim filepath As String = Path.Combine(Path.GetDirectoryName(System.Reflection.Assembly.GetExecutingAssembly().Location), "ReaktoroPropertyPackage.SetschenowCoefficients.txt")

        Using filestr As New FileStream(filepath, FileMode.OpenOrCreate)
            Using t As New IO.StreamReader(filestr)
                map = fh1.ReadStream(t)
            End Using
        End Using

        Maps.Clear()

        For Each item In map
            If Not Maps.ContainsKey(item.Name) Then Maps.Add(item.Name, item)
        Next

    End Sub

    Public Function GetValue(Name As String) As Double
        If Maps.ContainsKey(Name.ToLower()) Then
            Return Maps(Name.ToLower()).Coefficient
        Else
            Return 0.0
        End If
    End Function

End Class

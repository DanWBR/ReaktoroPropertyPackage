Imports System.Runtime.InteropServices

Module Module1

    <DllImport("kernel32.dll", SetLastError:=True)> Public Function AddDllDirectory(lpPathName As String) As Boolean

    End Function

    Public PythonPathSet As Boolean = False

End Module

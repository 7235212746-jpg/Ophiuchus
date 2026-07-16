using System;
using System.Diagnostics;
using System.IO;
using System.Windows.Forms;

internal static class Program
{
    [STAThread]
    private static void Main()
    {
        string baseDirectory = AppDomain.CurrentDomain.BaseDirectory;
        string launcher = Path.Combine(baseDirectory, "start_ophiuchus.bat");
        if (!File.Exists(launcher))
        {
            MessageBox.Show(
                "start_ophiuchus.bat was not found beside Ophiuchus.exe.",
                "Ophiuchus launcher",
                MessageBoxButtons.OK,
                MessageBoxIcon.Error);
            return;
        }

        try
        {
            string commandProcessor = Environment.GetEnvironmentVariable("COMSPEC") ?? "cmd.exe";
            ProcessStartInfo startInfo = new ProcessStartInfo
            {
                FileName = commandProcessor,
                Arguments = "/d /s /c \"\"" + launcher + "\"\"",
                WorkingDirectory = baseDirectory,
                UseShellExecute = false,
                CreateNoWindow = true,
                WindowStyle = ProcessWindowStyle.Hidden,
            };
            Process.Start(startInfo);
        }
        catch (Exception error)
        {
            MessageBox.Show(
                "Ophiuchus could not be started.\n\n" + error.Message,
                "Ophiuchus launcher",
                MessageBoxButtons.OK,
                MessageBoxIcon.Error);
        }
    }
}

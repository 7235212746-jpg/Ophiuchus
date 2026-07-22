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
        string portableApp = Path.Combine(baseDirectory, @"runtime\OphiuchusApp.exe");
        string launcher = Path.Combine(baseDirectory, "start_ophiuchus.bat");
        if (!File.Exists(portableApp) && !File.Exists(launcher))
        {
            MessageBox.Show(
                "Neither runtime\\OphiuchusApp.exe nor start_ophiuchus.bat was found beside Ophiuchus.exe.",
                "Ophiuchus launcher",
                MessageBoxButtons.OK,
                MessageBoxIcon.Error);
            return;
        }

        try
        {
            ProcessStartInfo startInfo;
            if (File.Exists(portableApp))
            {
                startInfo = new ProcessStartInfo
                {
                    FileName = portableApp,
                    WorkingDirectory = baseDirectory,
                    UseShellExecute = false,
                    CreateNoWindow = true,
                    WindowStyle = ProcessWindowStyle.Normal,
                };
            }
            else
            {
                string commandProcessor = Environment.GetEnvironmentVariable("COMSPEC") ?? "cmd.exe";
                startInfo = new ProcessStartInfo
                {
                    FileName = commandProcessor,
                    Arguments = "/d /s /c \"\"" + launcher + "\"\"",
                    WorkingDirectory = baseDirectory,
                    UseShellExecute = false,
                    CreateNoWindow = true,
                    WindowStyle = ProcessWindowStyle.Hidden,
                };
            }
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

#include "MyForm.h"
#include "MessageBOX.h"

using namespace System;
using namespace System::Windows::Forms;


[STAThread]
void Main()
{
	Application::EnableVisualStyles();
	Application::SetCompatibleTextRenderingDefault(false);

	Interface::MyForm form;
	Application::Run(%form);
}
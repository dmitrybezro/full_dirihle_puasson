#pragma once

namespace Interface {

	using namespace System;
	using namespace System::ComponentModel;
	using namespace System::Collections;
	using namespace System::Windows::Forms;
	using namespace System::Data;
	using namespace System::Drawing;

	/// <summary>
	/// Сводка для MessageBOX
	/// </summary>
	public ref class MessageBOX : public System::Windows::Forms::Form
	{
	public:
		MessageBOX(void)
		{
			InitializeComponent();
			//
			//TODO: добавьте код конструктора
			//
		}

	protected:
		/// <summary>
		/// Освободить все используемые ресурсы.
		/// </summary>
		~MessageBOX()
		{
			if (components)
			{
				delete components;
			}
		}
	public: System::Windows::Forms::TextBox^  textBox1;
	public: System::Windows::Forms::CheckBox^  checkBox1;
	public: System::Windows::Forms::Label^  label2;
	public: System::Windows::Forms::CheckBox^  checkBox2;
	public: System::Windows::Forms::Label^  label3;
	public: System::Windows::Forms::Label^  label4;
	public: System::Windows::Forms::Button^  button1;
	public:
	protected:

	private:
		/// <summary>
		/// Обязательная переменная конструктора.
		/// </summary>
		System::ComponentModel::Container ^components;

#pragma region Windows Form Designer generated code
		/// <summary>
		/// Требуемый метод для поддержки конструктора — не изменяйте 
		/// содержимое этого метода с помощью редактора кода.
		/// </summary>
		void InitializeComponent(void)
		{
			this->textBox1 = (gcnew System::Windows::Forms::TextBox());
			this->checkBox1 = (gcnew System::Windows::Forms::CheckBox());
			this->label2 = (gcnew System::Windows::Forms::Label());
			this->checkBox2 = (gcnew System::Windows::Forms::CheckBox());
			this->label3 = (gcnew System::Windows::Forms::Label());
			this->label4 = (gcnew System::Windows::Forms::Label());
			this->button1 = (gcnew System::Windows::Forms::Button());
			this->SuspendLayout();
			// 
			// textBox1
			// 
			this->textBox1->Location = System::Drawing::Point(213, 58);
			this->textBox1->Margin = System::Windows::Forms::Padding(4);
			this->textBox1->Name = L"textBox1";
			this->textBox1->Size = System::Drawing::Size(121, 26);
			this->textBox1->TabIndex = 0;
			this->textBox1->Text = L"8";
			// 
			// checkBox1
			// 
			this->checkBox1->AutoSize = true;
			this->checkBox1->Checked = true;
			this->checkBox1->CheckState = System::Windows::Forms::CheckState::Checked;
			this->checkBox1->Location = System::Drawing::Point(56, 117);
			this->checkBox1->Margin = System::Windows::Forms::Padding(4);
			this->checkBox1->Name = L"checkBox1";
			this->checkBox1->Size = System::Drawing::Size(264, 23);
			this->checkBox1->TabIndex = 1;
			this->checkBox1->Text = L"Установить оптимальное значение";
			this->checkBox1->UseVisualStyleBackColor = true;
			this->checkBox1->CheckedChanged += gcnew System::EventHandler(this, &MessageBOX::checkBox1_CheckedChanged);
			// 
			// label2
			// 
			this->label2->AutoSize = true;
			this->label2->Location = System::Drawing::Point(66, 94);
			this->label2->Margin = System::Windows::Forms::Padding(4, 0, 4, 0);
			this->label2->Name = L"label2";
			this->label2->Size = System::Drawing::Size(240, 19);
			this->label2->TabIndex = 3;
			this->label2->Text = L"Метод не имеет спец. параметров";
			// 
			// checkBox2
			// 
			this->checkBox2->AutoSize = true;
			this->checkBox2->Location = System::Drawing::Point(109, 147);
			this->checkBox2->Name = L"checkBox2";
			this->checkBox2->Size = System::Drawing::Size(176, 23);
			this->checkBox2->TabIndex = 4;
			this->checkBox2->Text = L"Сортировка значений";
			this->checkBox2->UseVisualStyleBackColor = true;
			// 
			// label3
			// 
			this->label3->AutoSize = true;
			this->label3->Location = System::Drawing::Point(96, 173);
			this->label3->Name = L"label3";
			this->label3->Size = System::Drawing::Size(210, 19);
			this->label3->TabIndex = 5;
			this->label3->Text = L"(Только для степеней двойки)";
			// 
			// label4
			// 
			this->label4->AutoSize = true;
			this->label4->Location = System::Drawing::Point(52, 58);
			this->label4->Name = L"label4";
			this->label4->Size = System::Drawing::Size(154, 19);
			this->label4->TabIndex = 6;
			this->label4->Text = L"Параметр Tau метода:";
			// 
			// button1
			// 
			this->button1->Location = System::Drawing::Point(168, 207);
			this->button1->Name = L"button1";
			this->button1->Size = System::Drawing::Size(55, 32);
			this->button1->TabIndex = 7;
			this->button1->Text = L"OK";
			this->button1->UseVisualStyleBackColor = true;
			this->button1->Click += gcnew System::EventHandler(this, &MessageBOX::button1_Click);
			// 
			// MessageBOX
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(9, 19);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->BackColor = System::Drawing::Color::FloralWhite;
			this->ClientSize = System::Drawing::Size(386, 251);
			this->Controls->Add(this->button1);
			this->Controls->Add(this->label4);
			this->Controls->Add(this->label3);
			this->Controls->Add(this->checkBox2);
			this->Controls->Add(this->label2);
			this->Controls->Add(this->checkBox1);
			this->Controls->Add(this->textBox1);
			this->Font = (gcnew System::Drawing::Font(L"Times New Roman", 12, System::Drawing::FontStyle::Regular, System::Drawing::GraphicsUnit::Point,
				static_cast<System::Byte>(204)));
			this->Margin = System::Windows::Forms::Padding(4);
			this->Name = L"MessageBOX";
			this->StartPosition = System::Windows::Forms::FormStartPosition::CenterParent;
			this->Text = L"Спец. параметры";
			this->ResumeLayout(false);
			this->PerformLayout();

		}
#pragma endregion
	private: System::Void button1_Click(System::Object^  sender, System::EventArgs^  e) {

		this->Hide();

	}
	private: System::Void checkBox1_CheckedChanged(System::Object^  sender, System::EventArgs^  e) {

		if (checkBox1->Checked)
			textBox1->ReadOnly = true;
		else textBox1->ReadOnly = false;

	}
};
}

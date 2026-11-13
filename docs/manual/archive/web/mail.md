# Mail

E-mail is the most common way to send messages to one or several addressees. It allows to the receiver to reply to the message or to forward it to others. It also allows to attach one or several files along. These messages can be stored and archived in the e-mail server or client.

## mail server

Service that runs in the background and manages user's emails account.

## mail client

Client application that fetches the client's e-mails from the server and sends outgoing e-mails to the server.

### Gmail IMAP client configuration

[This link](http://mail.google.com/support/bin/answer.py?answer=78799) explains how to configure your client to access your gmail with IMAP.

- Incoming Mail (IMAP) Server - requires SSL: imap.gmail.com
    -   Use SSL: Yes
    -   Port: 993
- Outgoing Mail (SMTP) Server - requires TLS: smtp.gmail.com (use authentication)
    -   Use Authentication: Yes
    -   Use STARTTLS: Yes (some clients call this SSL)
    -   Port: 465 or 587
- Account Name: your Gmail username (including @gmail.com)
- Email Address: your full Gmail email address (username@gmail.com)
- Password: your Gmail password

### Fenix IMAP client configuration

Para um cliente genérico configure o seguinte:

- Name: Coloque o seu nome como quer que apareça nos emails que envia.
- E-mail Address: Coloque o seu email. Por exemplo, antonio.silva@ist.utl.pt
- Incoming mail (IMAP): Preencha com mail.ist.utl.pt
    -   IMAP Use SSL: Selecione esta opção. Escolha o porto 995.
    -   Account Name: Coloque o login da sua conta de email. Por exemplo istxxxxx.
- Outgoing mail (SMTP): Preencha com mail.ist.utl.pt
    -   SMTP Use SSL: Selecione esta opção. Pode usar o porto 25 ou 465.
    -   SMTP authentication: Selecione esta opção. Selecione também a opção de usar os mesmos dados da conta de IMAP.

## External Links

- [Fenix webmail client configuration](https://ciist.ist.utl.pt/servicos/mail_config/index.php);
- [Gmail IMAP client configuration](http://mail.google.com/support/bin/answer.py?answer=78799).

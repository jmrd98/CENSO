import base64
import mimetypes
import urllib.parse

def file_to_base64_embed(file_path):
    """
    Converts a file to a Base64-encoded HTML <a> tag for embedding.

    Args:
        file_path (str): Path to the file to be encoded.

    Returns:
        str: JavaScript URI with document.write() to embed the Base64-encoded data.
    """
    try:
        # Determine the MIME type of the file
        mime_type, _ = mimetypes.guess_type(file_path)
        if mime_type not in ['application/pdf', 'image/png']:
            raise ValueError("Unsupported file type. Only PDF and PNG are supported.")
        
        # Read the file in binary mode
        with open(file_path, 'rb') as file:
            file_data = file.read()
        
        # Encode the binary data to Base64
        base64_data = base64.b64encode(file_data).decode('utf-8')
        
        # Split the base64 string into multiple lines
        line_length = 76
        base64_lines = [base64_data[i:i+line_length] for i in range(0, len(base64_data), line_length)]
        base64_string = '\n'.join(base64_lines)
        
        # Set the appropriate file extension for downloading
        extension = 'pdf' if mime_type == 'application/pdf' else 'png'
        
        # Create a complete HTML document with the Base64 data
        html_embed = f"""
        data:text/html,
        <!DOCTYPE html>
        <html>
        <head>
            <title>Base64 File Embed</title>
        </head>
        <body>
            <a href="data:{mime_type};base64,{base64_string}" download="file.{extension}">Download {extension.upper()}</a>
        </body>
        </html>
        """
        
        # Encode the HTML content for use in a URI
        encoded_html = urllib.parse.quote(html_embed)
        
        # Create the JavaScript URI
        js_uri = f"javascript:document.write(decodeURIComponent('{encoded_html}'))"
        
        return html_embed
    except Exception as e:
        return f"Error: {e}"

# Example usage
pdf_path = "nmrplot.pdf"
png_path = "nmrplot.png"

pdf_js_uri = file_to_base64_embed(pdf_path)
print(pdf_js_uri)

# png_js_uri = file_to_base64_embed(png_path)
# print(png_js_uri)
